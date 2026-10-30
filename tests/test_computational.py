"""Computational provenance contracts; all measurements here are synthetic."""
import csv
import json
import logging
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/lib"))
from computational import ComputationalRecorder, read_benchmark, render_report


def recorder(tmp):
    source = tmp / "Snakefile"
    source.write_text("# fixture")
    workflow = SimpleNamespace(dryrun=False, touch=False, global_resources={"_cores": 2, "mem_mb": 2048})
    result = ComputationalRecorder(workflow, tmp, source, "{}", ["sample"])
    result.start()
    return result


def start(rec, path, jobid=1):
    rec.job_started(SimpleNamespace(jobid=jobid, rule_name="fixture", wildcards={"isolate": "sample"},
        threads=2, resources={"mem_mb": 512}, benchmark=str(path) if path else None))


def benchmark(path, cpu="3.0", rss="12.0"):
    path.write_text("s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
                    f"2.0\t0:00:02\t{rss}\t20\t10\t11\t1.5\t0.5\t1.5\t{cpu}\n")


@pytest.mark.parametrize("value", ["NA", "-", "", "nan", "inf", "-1"])
def test_missing_metrics_are_not_zero(tmp_path, value):
    path = tmp_path / "benchmark.tsv"
    benchmark(path, value, value)
    row = read_benchmark(path)[0]
    assert row["cpu_seconds"] is None
    assert row["max_rss_mib"] is None
    assert row["io_read_mib"] == 1.5


def test_failed_retry_keeps_attempts_and_does_not_reuse_stale_benchmark(tmp_path):
    rec = recorder(tmp_path)
    path = tmp_path / "benchmark.tsv"
    benchmark(path, "999")  # Previous invocation must never count as new work.
    start(rec, path)
    rec.job_finished(1, "failed")
    start(rec, path)
    benchmark(path)
    rec.job_finished(1, "success")
    rec.finish("success")
    data = json.loads((rec.run_dir / "run.json").read_text())
    assert [j["attempt"] for j in data["jobs"]] == [1, 2]
    assert [j["status"] for j in data["jobs"]] == ["failed", "success"]
    assert data["jobs"][0]["measurements"] == []
    assert data["totals"]["measured_job_cpu_seconds_sum"] == 3
    assert data["totals"]["largest_job_peak_rss_mib"] == 12
    archive = rec.run_dir / data["jobs"][1]["benchmark_archive"]
    benchmark(path, "777")
    assert read_benchmark(archive)[0]["cpu_seconds"] == 3
    assert len((rec.run_dir / "events.jsonl").read_text().splitlines()) == 4
    rows = list(csv.DictReader((rec.run_dir / "jobs.tsv").open(), delimiter="\t"))
    assert rows[0]["cpu_seconds"] == "NA"
    assert rows[1]["cpu_seconds"] == "3.0"


def test_stale_or_absent_measurements_and_parallel_totals(tmp_path):
    rec = recorder(tmp_path)
    path = tmp_path / "old.tsv"
    benchmark(path)
    start(rec, path)
    rec.job_finished(1, "success")
    assert rec.data["jobs"][0]["measurement_status"] == "unavailable_no_fresh_benchmark"
    for jobid in (2, 3):
        path = tmp_path / f"{jobid}.tsv"
        start(rec, path, jobid)
        benchmark(path, "5", "100")
        rec.job_finished(jobid, "success")
    rec.finish("success")
    assert rec.data["totals"]["measured_job_cpu_seconds_sum"] == 10
    assert rec.data["totals"]["largest_job_peak_rss_mib"] == 100
    assert rec.data["elapsed_seconds"] < 2  # Never sum job wall times.


def test_quiet_completion_reconciles_dag_and_preserves_null_timing(tmp_path):
    rec = recorder(tmp_path)
    path = tmp_path / "job.tsv"
    start(rec, path)
    benchmark(path)
    rec.workflow.dag = SimpleNamespace(finished_jobs=[1], jobid=lambda job: job)
    rec.finish("success")
    job = rec.data["jobs"][0]
    assert job["status"] == "success"
    assert job["controller_elapsed_seconds"] is None
    assert job["measurements"][0]["wall_seconds"] == 2


def test_incomplete_exit_and_new_invocation_keep_history(tmp_path):
    rec = recorder(tmp_path)
    start(rec, None)
    rec.close()
    first = (rec.run_dir / "run.json").read_bytes()
    assert rec.data["status"] == "incomplete"
    assert rec.data["jobs"][0]["status"] == "incomplete"
    resumed = recorder(tmp_path)
    resumed.finish("no_work")
    assert (rec.run_dir / "run.json").read_bytes() == first
    assert resumed.data["totals"]["job_attempts"] == 0
    assert rec.run_dir != resumed.run_dir


def test_dryrun_and_touch_create_no_execution_records(tmp_path):
    for dryrun, touch in [(True, False), (False, True)]:
        workflow = SimpleNamespace(dryrun=False, touch=False, executor_plugin=None)
        rec = ComputationalRecorder(workflow, tmp_path, tmp_path / "absent", "{}", [])
        event = logging.makeLogRecord({"event": "workflow_started"})
        rec.emit(event)
        assert rec.data is None
        workflow.dryrun, workflow.touch, workflow.executor_plugin = dryrun, touch, object()
        rec.emit(logging.makeLogRecord({"msg": "Building DAG of jobs..."}))
        rec.close()
    assert not (tmp_path / "provenance").exists()


def test_report_escapes_labels_and_handles_empty_records(tmp_path):
    rec = recorder(tmp_path)
    start(rec, None)
    rec.data["jobs"][0]["rule"] = "<script>alert(1)</script>"
    rec.finish("failed")
    report = (rec.runtime_dir / "computational_report.html").read_text()
    assert "<script>" not in report
    assert "&lt;script&gt;" in report
    assert "unavailable" in report
    render_report(tmp_path / "empty")
