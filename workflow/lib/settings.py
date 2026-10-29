"""Resolve local run settings once; preserve explicit choices and resume identity.

This module selects existing presets and resource allocations. It does not infer
scientific thresholds, taxonomy, sequencing chemistry or new execution routes.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


POLICY_VERSION = 1
BASELINE = {
    "keep_percent": 95, "racon_rounds": 2, "busco_lineage": "bacteria_odb10",
    "medaka_model": "", "target_taxid": "", "allow_ont_only": True,
    "allow_illumina_only": True, "require_target_taxid_for_kraken": False,
    "kraken_decision_mode": "manual", "min_fasta_length": 1000,
}
MODULES = ("run_kraken", "run_busco", "run_prokka", "run_multiqc", "run_polypolish")


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False, default=str)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def is_auto(value):
    return isinstance(value, str) and value.strip().lower() == "auto"


def positive_int(value, name):
    if isinstance(value, bool) or not re.fullmatch(r"[0-9]+", str(value).strip()):
        raise ValueError(f"{name} must be a positive integer or auto, got {value!r}")
    result = int(value)
    if result < 1:
        raise ValueError(f"{name} must be positive")
    return result


def boolean(value, name):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        if value.strip().lower() in {"true", "yes", "y", "on", "1"}:
            return True
        if value.strip().lower() in {"false", "no", "n", "off", "0"}:
            return False
    raise ValueError(f"{name} must be auto, true/on or false/off, got {value!r}")


def _read(path):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return ""


def _cpuset_count(text):
    cpus = set()
    for part in text.split(","):
        if re.fullmatch(r"\d+(?:-\d+)?", part):
            bounds = list(map(int, part.split("-")))
            if bounds[-1] >= bounds[0] and bounds[-1] < 1_000_000:
                cpus.update(range(bounds[0], bounds[-1] + 1))
    return len(cpus) or None


def cgroup_limits(proc_root=Path("/proc")):
    """Read the process's v1/v2 cgroups and ancestor limits, including namespaces."""
    memberships = []
    for line in _read(Path(proc_root) / "self/cgroup").splitlines():
        fields = line.split(":", 2)
        if len(fields) == 3:
            memberships.append((set(fields[1].split(",")), fields[2]))
    limits = {"memory_bytes": [], "cpu_quotas": [], "cpuset_counts": []}
    for line in _read(Path(proc_root) / "self/mountinfo").splitlines():
        before, sep, after = line.partition(" - ")
        parts, info = before.split(), after.split()
        if not sep or len(parts) < 5 or len(info) < 3 or info[0] not in {"cgroup", "cgroup2"}:
            continue
        unescape = lambda s: re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), s)
        root, mount = unescape(parts[3]), Path(unescape(parts[4]))
        controllers = set(info[2].split(","))
        for member_controllers, member in memberships:
            v2 = info[0] == "cgroup2"
            if (v2 and member_controllers != {""}) or (not v2 and not controllers & member_controllers):
                continue
            if member == "/":
                relative = ""  # The cgroup namespace can expose its own root as /.
            elif root == "/":
                relative = member.lstrip("/")
            elif member == root or member.startswith(root.rstrip("/") + "/"):
                relative = member[len(root):].lstrip("/")
            else:
                continue
            if ".." in Path(relative).parts:
                continue
            current = mount / relative
            while True:
                if v2 or "memory" in controllers:
                    text = _read(current / ("memory.max" if v2 else "memory.limit_in_bytes"))
                    if text.isdigit() and 0 < int(text) < 2**60:
                        limits["memory_bytes"].append(int(text))
                if v2 or "cpu" in controllers:
                    quota = _read(current / "cpu.max").split() if v2 else [
                        _read(current / "cpu.cfs_quota_us"), _read(current / "cpu.cfs_period_us")]
                    if len(quota) == 2 and all(x.isdigit() and int(x) > 0 for x in quota):
                        limits["cpu_quotas"].append(int(quota[0]) / int(quota[1]))
                if v2 or "cpuset" in controllers:
                    count = _cpuset_count(_read(current / ("cpuset.cpus.effective" if v2 else "cpuset.cpus")))
                    if count:
                        limits["cpuset_counts"].append(count)
                if current == mount:
                    break
                current = current.parent
    return limits


def detect_capacity():
    logical = os.cpu_count() or 1
    try:
        affinity = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity = logical
    memory = None
    try:
        import psutil
        memory = int(psutil.virtual_memory().total)
    except (ImportError, OSError):
        try:
            memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (AttributeError, OSError, ValueError):
            pass
    limits = cgroup_limits()
    cpus = [logical, affinity, *limits["cpuset_counts"]]
    cpus += [max(1, int(quota)) for quota in limits["cpu_quotas"]]
    memories = ([memory] if memory and memory > 0 else []) + limits["memory_bytes"]
    return {"logical_cpus": logical, "affinity_cpus": affinity,
            "usable_cpus": max(1, min(cpus)), "physical_memory_bytes": memory,
            "memory_limit_bytes": min(memories) if memories else None,
            "cgroup_limits": limits}


def resource_value(resources, key):
    value = resources.get(key)
    return getattr(value, "value", value)


def kraken_memory_estimate(config, base):
    db = Path(base) / str(config.get("db_root", "db")) / str(config.get("kraken_db_name", "standard"))
    size = sum(p.stat().st_size for p in (db / f for f in ("hash.k2d", "opts.k2d", "taxo.k2d")) if p.is_file())
    return max(1024, int(size / 1048576 * 1.2) + 1024)


class SettingsResolution:
    def __init__(self, requested, *, cores=None, mem_mb=None, capacity=None):
        self.requested = json.loads(canonical(requested))
        self.capacity = capacity if capacity is not None else detect_capacity()
        mem_mb = positive_int(mem_mb, "Snakemake mem_mb") if mem_mb is not None else None
        self.scheduler = {"cores": cores, "mem_mb": mem_mb}
        base = Path(str(requested.get("base_dir", "."))).expanduser().absolute()
        sheet = requested.get("sample_sheet")
        sheet_path = Path(str(sheet)).absolute() if sheet else None
        context = {"base_dir": str(base), "working_directory": str(Path.cwd()),
                   "sample_sheet": str(sheet_path) if sheet_path else None,
                   "sample_sheet_sha256": hashlib.sha256(sheet_path.read_bytes()).hexdigest()
                   if sheet_path and sheet_path.is_file() else None}
        manifest = requested.get("input_manifest")
        if manifest:
            manifest_path = Path(str(manifest)).expanduser().absolute()
            context.update(input_manifest=str(manifest_path),
                           input_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest())
        self.identity = {"policy_version": POLICY_VERSION, "requested_config": self.requested,
                         "request_context": context}
        self.settings_id = digest(self.identity)
        self.path = base / "provenance/runtime/resolved_settings" / (self.settings_id + ".json")
        self.reused = self.path.exists()
        if self.reused:
            self.record = self._load()
        else:
            self.record = {"schema_version": 1, **deepcopy(self.identity), "settings_id": self.settings_id,
                           "created_at": datetime.now(timezone.utc).isoformat(),
                           "capacity_at_resolution": self.capacity, "scheduler_at_resolution": self.scheduler,
                           "resolved_config": deepcopy(self.requested), "decisions": {}, "sample_plans": {}}
            self._resolve()
        self.config = deepcopy(self.record["resolved_config"])
        if (self.reused and self.config["run_kraken"] and self.config["kraken_decision_mode"] != "off"
                and is_auto(self.requested.get("kraken_mem_mb", "auto"))
                and kraken_memory_estimate(self.config, base) != self.config["kraken_mem_mb"]):
            raise ValueError("Kraken database size changes the saved Auto memory estimate. Use a new settings_revision to resolve its current requirement; keep the old record for provenance.")
        ceilings = [self.config["total_mem_mb"]]
        if mem_mb is not None:
            ceilings.append(positive_int(mem_mb, "Snakemake mem_mb"))
        if self.capacity["memory_limit_bytes"]:
            fraction = 0.8 if is_auto(self.requested.get("total_mem_mb", "auto")) else 1
            ceilings.append(int(self.capacity["memory_limit_bytes"] * fraction) // 1_000_000)
        self.memory_budget = min(ceilings)
        self.cpu_budget = min(self.capacity["usable_cpus"], positive_int(cores, "Snakemake cores")) if cores else self.capacity["usable_cpus"]

    def _load(self):
        try:
            record = json.loads(self.path.read_text())
            checksum = record.pop("record_sha256")
            if checksum != digest(record) or any(record[k] != v for k, v in self.identity.items()):
                raise ValueError("identity or checksum mismatch")
            if record["schema_version"] != 1 or record["settings_id"] != self.settings_id:
                raise ValueError("unsupported record")
            for key in ("threads", "total_mem_mb", "assembly_mem_mb", "spades_memory_gb",
                        "polishing_mem_mb", "kraken_mem_mb", "busco_mem_mb", "prokka_mem_mb"):
                positive_int(record["resolved_config"][key], key)
            return record
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise ValueError(f"Cannot reuse resolved settings {self.path}: {exc}. Preserve this record and review it; do not silently re-detect a resumed run.") from exc

    def _choose(self, key, automatic, reason):
        requested = self.requested.get(key, "auto")
        auto = is_auto(requested)
        value = automatic if auto else requested
        self.record["resolved_config"][key] = value
        self.record["decisions"][key] = {"requested": requested, "value": value,
            "source": "auto" if auto else "manual", "reason": reason if auto else "Explicit configuration override"}
        return value

    def _integer(self, key, automatic, reason):
        value = positive_int(self._choose(key, automatic, reason), key)
        self.record["resolved_config"][key] = value
        self.record["decisions"][key]["value"] = value
        return value

    def _resolve(self):
        cfg = self.record["resolved_config"]
        for key, value in BASELINE.items():
            if key not in cfg:
                cfg[key] = value
                self.record["decisions"][key] = {"requested": None, "value": value,
                    "source": "preset", "reason": "Existing v9 processing default"}
            else:
                self.record["decisions"][key] = {"requested": cfg[key], "value": cfg[key],
                    "source": "manual", "reason": "Explicit configuration override"}
        cpus = self.capacity["usable_cpus"]
        if self.scheduler["cores"]:
            cpus = min(cpus, positive_int(self.scheduler["cores"], "Snakemake cores"))
        self._integer("threads", cpus, "Available CPU affinity/cgroup capacity bounded by Snakemake cores")
        memory = self.capacity["memory_limit_bytes"]
        manual_budget = self.scheduler["mem_mb"]
        if memory:
            automatic_budget = max(1, int(memory * 0.8) // 1_000_000)
            if manual_budget is not None:
                automatic_budget = min(automatic_budget, positive_int(manual_budget, "Snakemake mem_mb"))
        elif manual_budget is not None:
            automatic_budget = positive_int(manual_budget, "Snakemake mem_mb")
        elif is_auto(self.requested.get("total_mem_mb", "auto")):
            raise ValueError("Cannot detect available memory; provide total_mem_mb or Snakemake --resources mem_mb=...")
        else:
            automatic_budget = 1  # Unused: the explicit value takes precedence.
        budget = self._integer("total_mem_mb", automatic_budget,
            "80% of the physical/cgroup memory limit, bounded by Snakemake mem_mb" if memory else "Explicit Snakemake memory budget; host memory unavailable")
        if memory and budget > memory // 1_000_000:
            raise ValueError("total_mem_mb exceeds the detected host/container memory limit")
        usable_budget = min(budget, manual_budget) if manual_budget is not None else budget
        spades = self.requested.get("spades_memory_gb", "auto")
        required = 1024 if is_auto(spades) else positive_int(spades, "spades_memory_gb") * 1024
        assembly = self._integer("assembly_mem_mb", max(required, min(8192, usable_budget)),
            "Existing 8192 MB reservation bounded by total memory; respect explicit SPAdes memory")
        spades = self._integer("spades_memory_gb", max(1, assembly // 1024), "Fit within the resolved assembly reservation")
        if spades * 1024 > assembly:
            raise ValueError("spades_memory_gb must fit within assembly_mem_mb")
        for key in ("polishing_mem_mb", "busco_mem_mb", "prokka_mem_mb"):
            self._integer(key, min(4096, usable_budget), "Existing 4096 MB reservation bounded by total memory")
        self._integer("kraken_mem_mb", kraken_memory_estimate(cfg, self.identity["request_context"]["base_dir"]),
                      "Existing database-size estimate; never shrink it to conceal a capacity shortfall")
        for key in MODULES:
            enabled = not (key == "run_kraken" and str(cfg["kraken_decision_mode"]).strip().lower() == "off")
            value = boolean(self._choose(key, enabled,
                "Existing module default; input eligibility is recorded per sample" if enabled else "Kraken policy is explicitly off"), key)
            cfg[key] = self.record["decisions"][key]["value"] = value
        value = str(self._choose("polish_mode", "tripolish", "Existing Tripolish preset; only input-compatible stages are planned")).strip().lower()
        if value not in {"none", "racon", "medaka", "tripolish"}:
            raise ValueError("polish_mode must be auto, none, racon, medaka or tripolish")
        cfg["polish_mode"] = self.record["decisions"]["polish_mode"]["value"] = value
        mode = str(cfg.get("run_mode", "auto")).strip().lower()
        if mode not in {"auto", "hybrid", "illumina_only", "ont_only"}:
            raise ValueError("run_mode must be auto, hybrid, illumina_only or ont_only")
        cfg["run_mode"] = mode

    def set_sample_plans(self, plans):
        if self.reused and self.record["sample_plans"] != plans:
            raise ValueError("Inputs or sample settings differ from the saved resolved plan. Use a new draft or change settings_revision explicitly to resolve a new plan; keep the old record for provenance.")
        self.record["sample_plans"] = deepcopy(plans)

    def persist(self):
        """Publish a complete record atomically, without replacing any prior plan."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        record = {**self.record, "record_sha256": digest(self.record)}
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=self.path.parent, delete=False) as tmp:
            temporary = Path(tmp.name)
            json.dump(record, tmp, indent=2, sort_keys=True, allow_nan=False)
            tmp.write("\n")
            tmp.flush()
            os.fsync(tmp.fileno())
        try:
            try:
                os.link(temporary, self.path)
            except FileExistsError:
                existing = self._load()
                if any(existing[k] != self.record[k] for k in ("resolved_config", "decisions", "sample_plans")):
                    raise ValueError("Another controller resolved different settings for this request; retry to reuse its saved plan")
                self.record = existing
            self.reused = True
        finally:
            temporary.unlink(missing_ok=True)

    def summary(self):
        return (f"[HybAs settings] {'Reusing' if self.reused else 'Previewing'} {self.settings_id[:12]}: "
                f"threads={self.config['threads']}, scheduler cores={self.cpu_budget}, "
                f"memory budget={self.memory_budget} MB, polish_mode={self.config['polish_mode']}. "
                f"Resolved settings: {self.path}")
