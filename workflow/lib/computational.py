"""Persist computational provenance from Snakemake's controller events.

Job resource measurements come from native benchmarks, not from resource requests.
This module never changes scheduling or analysis outputs.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import logging
import math
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import uuid


METRICS = {
    "s": "wall_seconds", "cpu_time": "cpu_seconds", "max_rss": "max_rss_mib",
    "max_vms": "max_vms_mib", "max_uss": "max_uss_mib", "max_pss": "max_pss_mib",
    # Snakemake 9.27.1 divides the psutil byte counters by 1024**2.
    "io_in": "io_read_mib", "io_out": "io_write_mib", "mean_load": "mean_cpu_load",
    "input_size_mb": "input_size_mib",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_json(path, value):
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) and result >= 0 else None
    except (TypeError, ValueError):
        return None


def fingerprint(path):
    try:
        stat = Path(path).stat()
        return stat.st_mtime_ns, stat.st_size
    except (OSError, TypeError):
        return None


def read_benchmark(path):
    """Keep unavailable measurements null, including on very short jobs."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "s" not in reader.fieldnames:
            raise ValueError("Benchmark has no wall-time column")
        rows = [{metric: number(row.get(column)) for column, metric in METRICS.items()}
                for row in reader]
    if not rows:
        raise ValueError("Benchmark contains no measurements")
    return rows


def machine_info(base_dir):
    info = dict(system=platform.system(), release=platform.release(),
                architecture=platform.machine(), python=platform.python_version(),
                logical_cpus=os.cpu_count(), physical_cpus=None, ram_total_bytes=None,
                cpu_affinity_count=None, disk_free_bytes=shutil.disk_usage(base_dir).free)
    try:
        from snakemake import __version__
        info["snakemake"] = __version__
    except ImportError:
        info["snakemake"] = None
    try:
        import psutil
        info["physical_cpus"] = psutil.cpu_count(logical=False)
        info["ram_total_bytes"] = psutil.virtual_memory().total
        info["cpu_affinity_count"] = len(psutil.Process().cpu_affinity())
    except Exception as error:
        # /proc and CPU affinity are not available on every execution host.
        info["hardware_collection_note"] = str(error)
    return info


def plain_resources(resources):
    return {str(key): getattr(value, "value", value) for key, value in resources.items()
            if isinstance(getattr(value, "value", value), (str, int, float, bool, type(None)))
            and not (isinstance(getattr(value, "value", value), float)
                     and not math.isfinite(getattr(value, "value", value)))}


def totals(jobs):
    measurements = [row for job in jobs for row in job.get("measurements", [])]
    cpu = [row["cpu_seconds"] for row in measurements if row["cpu_seconds"] is not None]
    rss = [row["max_rss_mib"] for row in measurements if row["max_rss_mib"] is not None]
    return {
        "job_attempts": len(jobs),
        "attempts_expecting_benchmarks": sum(bool(job.get("benchmark_path")) for job in jobs),
        "successful_attempts": sum(job["status"] == "success" for job in jobs),
        "failed_attempts": sum(job["status"] == "failed" for job in jobs),
        "attempts_with_benchmarks": sum(bool(job.get("measurements")) for job in jobs),
        "attempts_with_cpu_measurements": sum(any(row["cpu_seconds"] is not None
            for row in job.get("measurements", [])) for job in jobs),
        "measured_job_cpu_seconds_sum": sum(cpu) if cpu else None,
        "largest_job_peak_rss_mib": max(rss) if rss else None,
    }


def render_report(runtime_dir):
    """Standalone report, embedded by Snakemake; data remains in immutable run dirs."""
    runs = []
    for path in sorted((Path(runtime_dir) / "computational_runs").glob("*/run.json"), reverse=True):
        try:
            runs.append(json.loads(path.read_text()))
        except (OSError, ValueError):
            continue

    def cell(value):
        if value is None:
            return "unavailable"
        return html.escape(f"{value:.3f}" if isinstance(value, float) else str(value))

    def table(headers, rows):
        return "<table><thead><tr>" + "".join(f"<th>{cell(h)}</th>" for h in headers) + \
            "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{cell(v)}</td>" for v in row) +
            "</tr>" for row in rows) + "</tbody></table>"

    content = ["<!doctype html><html lang='en'><meta charset='utf-8'>",
        "<title>HybAs computational record</title><style>body{font:15px system-ui;max-width:1300px;margin:32px auto;padding:0 20px;color:#173047}table{border-collapse:collapse;width:100%;margin:20px 0}th,td{padding:9px;text-align:left;border-bottom:1px solid #dbe4eb}th{background:#edf4f8}details{margin:20px 0;overflow:auto}small{color:#526575}</style>",
        "<h1>Computational record</h1><p>One record per execution invocation. Earlier runs and retries are retained.</p>",
        "<p>Elapsed time includes planning and scheduling. CPU time and memory come from sampled native job benchmarks. The largest job peak is <strong>not</strong> simultaneous workflow peak RAM. Missing values are unavailable, not zero. Rule success includes conditional rules that may skip their tool.</p>"]
    content.append(table(["Run (UTC)", "Status", "Elapsed seconds", "Attempts", "Measured CPU seconds", "Largest job RSS (MiB)"],
        [[r["started_at"], r["status"], r.get("elapsed_seconds"), r["totals"]["job_attempts"],
          r["totals"]["measured_job_cpu_seconds_sum"], r["totals"]["largest_job_peak_rss_mib"]] for r in runs]))
    for run in runs:
        content.append("<details><summary>" + cell(run["run_id"]) + " — " + cell(run["status"]) + "</summary>")
        if run.get("settings_decisions"):
            content.append("<details><summary>Resolved run settings</summary>")
            content.append(table(["Setting", "Requested", "Selected", "Source", "Reason"],
                [[name, d["requested"], d["value"], d["source"], d["reason"]]
                 for name, d in sorted(run["settings_decisions"].items())]))
            content.append(table(["Sample", "Input mode", "Mode reason", "Planned polishing stages"],
                [[sample, plan["mode"], plan["mode_reason"], ", ".join(plan["planned_polishing_stages"]) or "none"]
                 for sample, plan in sorted(run.get("sample_plans", {}).items())]))
            content.append("<small>Settings snapshot: " + cell(run["settings_id"]) + "</small></details>")
        content.append(table(["Rule", "Sample", "Attempt", "Status", "Threads requested", "Wall s", "CPU s", "Peak RSS MiB", "Read MiB", "Write MiB", "Measurement status"],
            [[j["rule"], j["wildcards"].get("isolate", "global"), j["attempt"], j["status"], j["threads_requested"],
              m.get("wall_seconds"), m.get("cpu_seconds"), m.get("max_rss_mib"), m.get("io_read_mib"),
              m.get("io_write_mib"), j.get("measurement_status", "pending")]
             for j in run["jobs"] for m in (j.get("measurements") or [{}])]))
        content.append("<small>Controller machine and runner: " + cell(json.dumps(run["machine"], sort_keys=True)) + "</small></details>")
    content.append("</html>")
    atomic_text(Path(runtime_dir) / "computational_report.html", "\n".join(content))


class ComputationalRecorder(logging.Handler):
    """Synchronous event capture: persist starts and endings before possible retry."""

    def __init__(self, workflow, base_dir, source, config_json, samples, settings=None):
        super().__init__()
        self.workflow = workflow
        self.base_dir = Path(base_dir)
        self.runtime_dir = self.base_dir / "provenance/runtime"
        self.source = Path(source)
        self.config_hash = hashlib.sha256(config_json.encode()).hexdigest()
        self.samples = list(samples)
        self.settings = settings
        self.owner_pid = os.getpid()
        self.data = None
        self.active = {}
        self.warned = False
        self.finalized = False
        self.waiting_for_executor = False

    def emit(self, record):
        if os.getpid() != self.owner_pid:
            return
        # ExecutionSettings does not exist while parsing the Snakefile. Check
        # the mandatory freeze hook once Snakemake has supplied those settings,
        # before DAG execution. This error must not become a metrics warning.
        if (self.settings is not None and str(getattr(record, "event", "")) == "workflow_started"
                and self.workflow.execution_settings.no_hooks):
            message = "HybAs requires start hooks to freeze resolved settings; omit --no-hooks"
            print(message, file=sys.stderr)  # Snakemake --quiet all also hides its error renderer.
            raise ValueError(message)
        try:
            event = str(getattr(record, "event", ""))
            if event == "workflow_started":
                # Snakemake emits this BEFORE selecting its executor. At that
                # instant workflow.dryrun and workflow.touch are both false.
                self.waiting_for_executor = True
                return
            if self.waiting_for_executor and self.workflow.executor_plugin is not None:
                self.waiting_for_executor = False
                if not self.workflow.dryrun and not self.workflow.touch:
                    self.start()
            if self.data is None or self.finalized:
                return
            if event == "job_info":
                self.job_started(record)
            elif event == "job_finished":
                self.job_finished(record.job_id, "success")
            elif event == "job_error":
                self.job_finished(record.jobid, "failed")
            elif event == "group_error":
                for job in record.job_error_info:
                    self.job_finished(job["jobid"], "failed")
            elif record.getMessage().startswith("Nothing to be done"):
                self.finish("no_work")
        except Exception as error:
            self.warn(error)

    def warn(self, error):
        if not self.warned:
            print(f"WARNING: HybAs computational recording is incomplete: {error}", file=sys.stderr)
            self.warned = True
        if self.data is not None:
            self.data.setdefault("collection_errors", []).append(str(error))
            try:
                self.save()
            except Exception:
                pass  # The warning remains visible even if storage itself failed.

    def start(self):
        if self.data is not None:
            return
        self.start_clock = time.monotonic()
        run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + "-" + uuid.uuid4().hex[:8]
        self.run_dir = self.runtime_dir / "computational_runs" / run_id
        self.run_dir.mkdir(parents=True)
        self.data = dict(schema_version=1, run_id=run_id, started_at=utc_now(), ended_at=None,
            status="running", elapsed_seconds=None, samples=self.samples, jobs=[],
            config_sha256=self.config_hash, source_sha256=hashlib.sha256(self.source.read_bytes()).hexdigest(),
            machine=machine_info(self.base_dir), requested_resources=plain_resources(self.workflow.global_resources),
            measurement_method="Snakemake native sampled job benchmarks", collection_errors=[])
        if self.settings is not None:
            self.data.update(settings_id=self.settings.settings_id,
                settings_snapshot=str(self.settings.path),
                settings_decisions=self.settings.record["decisions"],
                sample_plans=self.settings.record["sample_plans"])
            write_json(self.run_dir / "resolved_settings.json", self.settings.record)
        self.save()

    def event(self, kind, job):
        with (self.run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(dict(event=kind, timestamp=utc_now(), job=job), allow_nan=False) + "\n")
            handle.flush()

    def job_started(self, record):
        jobid = str(record.jobid)
        if jobid in self.active:
            return  # Duplicate group announcements are not extra attempts.
        attempt = 1 + sum(j["job_id"] == jobid for j in self.data["jobs"])
        job = dict(job_id=jobid, rule=record.rule_name, wildcards=dict(record.wildcards), attempt=attempt,
            status="running", started_at=utc_now(), ended_at=None, controller_elapsed_seconds=None,
            threads_requested=record.threads, resources_requested=plain_resources(record.resources),
            benchmark_path=str(record.benchmark) if record.benchmark else None, measurements=[])
        self.data["jobs"].append(job)
        self.active[jobid] = (job, time.monotonic(), fingerprint(job["benchmark_path"]))
        self.event("job_started", job)
        self.save()

    def job_finished(self, jobid, status, inferred=False):
        active = self.active.pop(str(jobid), None)
        if active is None:
            return
        job, clock, before = active
        job.update(status=status, ended_at=utc_now(),
                   completion_time_source="workflow_finalization" if inferred else "controller_event",
                   controller_elapsed_seconds=None if inferred else time.monotonic() - clock)
        path = job["benchmark_path"]
        if status != "success":
            job["measurement_status"] = "unavailable_failed_or_interrupted"
        elif not path:
            job["measurement_status"] = "not_benchmarked"
        elif fingerprint(path) is None or fingerprint(path) == before:
            job["measurement_status"] = "unavailable_no_fresh_benchmark"
        else:
            try:
                job["measurements"] = read_benchmark(path)
                relative = Path("benchmarks") / f"job-{job['job_id']}-attempt-{job['attempt']}.tsv"
                target = self.run_dir / relative
                target.parent.mkdir(exist_ok=True)
                shutil.copyfile(path, target)
                job["benchmark_archive"] = str(relative)
                job["measurement_status"] = "sampled" if all(
                    row["cpu_seconds"] is not None and row["max_rss_mib"] is not None
                    for row in job["measurements"]) else "partial"
            except (OSError, ValueError) as error:
                job["measurement_status"] = "unavailable_invalid_benchmark"
                job["measurement_error"] = str(error)
        self.event("job_finished", job)
        self.save()

    def save(self):
        self.data["totals"] = totals(self.data["jobs"])
        write_json(self.run_dir / "run.json", self.data)

    def finish(self, status):
        with self.lock:
            if self.data is None or self.finalized:
                return
            try:
                # --quiet progress suppresses JOB_FINISHED events in Snakemake.
                # Reconcile from the DAG, without inventing completion timestamps.
                dag = getattr(self.workflow, "dag", None)
                completed = {str(dag.jobid(job)) for job in dag.finished_jobs} if dag is not None else set()
                for jobid in list(self.active):
                    self.job_finished(jobid, "success" if jobid in completed else "incomplete", inferred=True)
                self.data.update(status=status, ended_at=utc_now(), elapsed_seconds=time.monotonic() - self.start_clock)
                self.save()
                fields = ["run_id", "job_id", "rule", "sample", "attempt", "status", "threads_requested",
                          "controller_elapsed_seconds", "measurement_status", *METRICS.values()]
                with (self.run_dir / "jobs.tsv").open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fields, delimiter="\t")
                    writer.writeheader()
                    for job in self.data["jobs"]:
                        for measurement in job.get("measurements") or [{}]:
                            row = {key: job.get(key) for key in fields if key not in METRICS.values()}
                            row.update(run_id=self.data["run_id"], sample=job["wildcards"].get("isolate"), **measurement)
                            writer.writerow({key: "NA" if row.get(key) is None else row[key] for key in fields})
                self.finalized = True
                render_report(self.runtime_dir)
            except Exception as error:
                self.warn(error)

    def close(self):
        # Snakemake shuts down handlers on exceptional exits too. A hard kill
        # cannot finalize; its already-written record stays explicitly running.
        if os.getpid() == self.owner_pid:
            self.finish("incomplete")
        super().close()


def install_recorder(workflow, base_dir, source, config_json, samples, settings=None):
    from snakemake.logging import logger
    workflow.output_settings.benchmark_extended = True
    if not workflow.is_main_process:
        return None
    recorder = ComputationalRecorder(workflow, base_dir, source, config_json, samples, settings=settings)
    logger.addHandler(recorder)
    return recorder
