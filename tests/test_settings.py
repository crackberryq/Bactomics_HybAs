"""Configuration/host fixtures only; no scientific tools or data."""
from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/lib"))
from settings import SettingsResolution, cgroup_limits, digest


def capacity(cpus=12, memory_mb=32000):
    return dict(usable_cpus=cpus, memory_limit_bytes=memory_mb * 1_000_000 if memory_mb else None)


def resolve(tmp, values=None, **kwargs):
    return SettingsResolution({"base_dir": str(tmp), **(values or {})},
                              capacity=kwargs.pop("capacity", capacity()), **kwargs)


def test_auto_uses_capacity_and_scheduler_and_preserves_existing_preset(tmp_path):
    s = resolve(tmp_path, cores=4, mem_mb=12000)
    assert s.config["threads"] == 4
    assert s.config["total_mem_mb"] == 12000
    assert s.config["assembly_mem_mb"] == 8192
    assert s.config["spades_memory_gb"] == 8
    assert s.config["polishing_mem_mb"] == 4096
    assert s.config["polish_mode"] == "tripolish"
    assert s.config["run_busco"] is True  # Missing databases must not silently disable it.
    assert s.record["decisions"]["threads"]["source"] == "auto"
    assert not s.path.exists()  # Preview is read-only.


def test_manifest_content_is_part_of_resume_identity(tmp_path):
    manifest = tmp_path / 'inputs.yaml'
    manifest.write_text('fixture: first\n')
    first = resolve(tmp_path, {'input_manifest': str(manifest)})
    first.persist()
    assert resolve(tmp_path, {'input_manifest': str(manifest)}).reused
    manifest.write_text('fixture: changed\n')
    changed = resolve(tmp_path, {'input_manifest': str(manifest)})
    assert changed.settings_id != first.settings_id
    assert not changed.reused and first.path.exists()


def test_small_machine_and_manual_spades_constraint(tmp_path):
    s = resolve(tmp_path, capacity=capacity(2, 4000))
    assert s.config["total_mem_mb"] == 3200
    assert s.config["assembly_mem_mb"] == 3200
    assert s.config["spades_memory_gb"] == 3
    assert s.config["polishing_mem_mb"] == 3200
    manual = resolve(tmp_path, {"spades_memory_gb": 12})
    assert manual.config["assembly_mem_mb"] == 12 * 1024
    assert manual.config["spades_memory_gb"] == 12


def test_explicit_values_win_without_silent_clipping(tmp_path):
    s = resolve(tmp_path, {"threads": 16, "total_mem_mb": 16000, "assembly_mem_mb": 10000,
        "spades_memory_gb": 7, "polishing_mem_mb": 2000, "polish_mode": "none",
        "run_busco": "off", "run_prokka": False, "keep_percent": 88}, cores=2, mem_mb=8000)
    assert s.config["threads"] == 16
    assert s.cpu_budget == 2
    assert s.config["total_mem_mb"] == 16000
    assert s.memory_budget == 8000
    assert s.config["assembly_mem_mb"] == 10000  # Workflow preflight rejects this reservation.
    assert s.config["spades_memory_gb"] == 7
    assert s.config["polishing_mem_mb"] == 2000
    assert s.config["polish_mode"] == "none"
    assert s.config["run_busco"] is False
    assert s.config["keep_percent"] == 88
    assert s.record["decisions"]["threads"]["source"] == "manual"


@pytest.mark.parametrize("key,value", [("threads", True), ("threads", 0), ("threads", 1.5),
    ("total_mem_mb", -1), ("total_mem_mb", None), ("run_busco", "maybe"),
    ("run_prokka", 1), ("polish_mode", "guess"), ("run_mode", "unsupported_future_mode")])
def test_invalid_explicit_values_do_not_become_auto(tmp_path, key, value):
    with pytest.raises(ValueError):
        resolve(tmp_path, {key: value})


def test_manual_memory_cannot_exceed_machine_or_conflict_with_spades(tmp_path):
    with pytest.raises(ValueError, match="detected host/container"):
        resolve(tmp_path, {"total_mem_mb": 33000})
    with pytest.raises(ValueError, match="fit within"):
        resolve(tmp_path, {"assembly_mem_mb": 1024, "spades_memory_gb": 2})


def test_kraken_auto_flag_respects_explicit_off_and_keeps_explicit_mode(tmp_path):
    s = resolve(tmp_path, {"run_kraken": "auto", "kraken_decision_mode": "off"})
    assert s.config["run_kraken"] is False
    s = resolve(tmp_path, {"run_kraken": False, "kraken_decision_mode": "aware"})
    assert s.config["run_kraken"] is False  # Existing workflow contradiction check remains effective.
    assert s.config["kraken_decision_mode"] == "aware"
    s = resolve(tmp_path, {"kraken_decision_mode": "auto"})
    assert s.config["kraken_decision_mode"] == "auto"  # Never reinterpret the existing Kraken algorithm.


def test_database_estimate_not_capped_to_host_budget(tmp_path):
    db = tmp_path / "db/standard"
    db.mkdir(parents=True)
    with (db / "hash.k2d").open("wb") as handle:
        handle.truncate(10 * 1024**3)  # Sparse fixture, no database data.
    s = resolve(tmp_path, mem_mb=4000)
    assert s.config["kraken_mem_mb"] > s.memory_budget


def test_resume_refuses_a_stale_auto_database_memory_estimate(tmp_path):
    first = resolve(tmp_path)
    first.persist()
    original = first.path.read_bytes()
    db = tmp_path / "db/standard"
    db.mkdir(parents=True)
    with (db / "hash.k2d").open("wb") as handle:
        handle.truncate(100 * 1024**2)
    with pytest.raises(ValueError, match="database size changes"):
        resolve(tmp_path)
    fresh = resolve(tmp_path, {"settings_revision": "updated-database"})
    assert fresh.config["kraken_mem_mb"] > first.config["kraken_mem_mb"]
    assert first.path.read_bytes() == original


def test_unknown_memory_requires_a_budget(tmp_path):
    with pytest.raises(ValueError, match="Cannot detect"):
        resolve(tmp_path, capacity=capacity(memory_mb=0))
    s = resolve(tmp_path, mem_mb=6000, capacity=capacity(memory_mb=0))
    assert s.memory_budget == 6000
    s = resolve(tmp_path, {"total_mem_mb": 7000}, capacity=capacity(memory_mb=0))
    assert s.memory_budget == 7000


def test_resume_reuses_choices_and_rechecks_current_ceiling(tmp_path):
    first = resolve(tmp_path, cores=6)
    first.set_sample_plans({"s": {"mode": "illumina_only"}})
    first.persist()
    before = first.path.read_bytes(), first.path.stat().st_mtime_ns
    again = resolve(tmp_path, cores=2, capacity=capacity(2, 16000))
    assert again.reused
    assert again.config == first.config
    assert again.config["threads"] == 6
    assert again.cpu_budget == 2
    assert again.memory_budget == 12800
    again.set_sample_plans({"s": {"mode": "illumina_only"}})
    again.persist()
    assert (again.path.read_bytes(), again.path.stat().st_mtime_ns) == before
    with pytest.raises(ValueError, match="Inputs or sample settings differ"):
        again.set_sample_plans({"s": {"mode": "hybrid"}})


def test_explicit_edit_revision_or_sample_sheet_change_gets_new_identity(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("isolate\ns1\n")
    values = {"sample_sheet": str(sheet)}
    first = resolve(tmp_path, values)
    first.persist()
    original = first.path.read_bytes()
    assert resolve(tmp_path, {**values, "threads": 3}).settings_id != first.settings_id
    assert resolve(tmp_path, {**values, "settings_revision": "new-host"}).settings_id != first.settings_id
    sheet.write_text("isolate\ns1\ns2\n")
    assert resolve(tmp_path, values).settings_id != first.settings_id
    assert first.path.read_bytes() == original


@pytest.mark.parametrize("broken", ["{", "{}", "tampered"])
def test_corrupt_snapshot_blocks_re_detection(tmp_path, broken):
    s = resolve(tmp_path)
    s.persist()
    if broken == "tampered":
        data = json.loads(s.path.read_text())
        data["resolved_config"]["threads"] += 1
        broken = json.dumps(data)
    s.path.write_text(broken)
    with pytest.raises(ValueError, match="Cannot reuse resolved settings"):
        resolve(tmp_path)


def test_concurrent_different_resolution_cannot_replace_first_snapshot(tmp_path):
    first = resolve(tmp_path, cores=2)
    second = resolve(tmp_path, cores=4)
    first.persist()
    original = first.path.read_bytes()
    with pytest.raises(ValueError, match="Another controller"):
        second.persist()
    assert first.path.read_bytes() == original
    assert len(list(first.path.parent.iterdir())) == 1


def test_record_checksum_and_original_request_are_saved(tmp_path):
    request = {"threads": "auto", "run_prokka": False}
    original = deepcopy(request)
    s = resolve(tmp_path, request)
    s.persist()
    assert request == original
    record = json.loads(s.path.read_text())
    checksum = record.pop("record_sha256")
    assert digest(record) == checksum
    assert record["requested_config"]["threads"] == "auto"
    assert record["resolved_config"]["threads"] == 12


def proc_fixture(tmp, lines, mounts):
    proc = tmp / "proc/self"
    proc.mkdir(parents=True)
    (proc / "cgroup").write_text(lines)
    (proc / "mountinfo").write_text(mounts)
    return proc.parent


def test_cgroup_v2_inherits_parent_limits(tmp_path):
    mount = tmp_path / "cgroup"
    child = mount / "parent/child"
    child.mkdir(parents=True)
    (child / "memory.max").write_text("max")
    (child.parent / "memory.max").write_text("4000000000")
    (mount / "memory.max").write_text("8000000000")
    (child / "cpu.max").write_text("max 100000")
    (child.parent / "cpu.max").write_text("250000 100000")
    (child / "cpuset.cpus.effective").write_text("2-4,7")
    proc = proc_fixture(tmp_path, "0::/parent/child\n", f"1 0 0:1 / {mount} rw - cgroup2 cgroup rw\n")
    result = cgroup_limits(proc)
    assert min(result["memory_bytes"]) == 4000000000
    assert result["cpu_quotas"] == [2.5]
    assert result["cpuset_counts"] == [4]


def test_cgroup_v1_mount_root_and_namespace(tmp_path):
    memory, cpu = tmp_path / "memory", tmp_path / "cpu"
    memory.mkdir()
    cpu.mkdir()
    (memory / "memory.limit_in_bytes").write_text("2000000000")
    (cpu / "cpu.cfs_quota_us").write_text("100000")
    (cpu / "cpu.cfs_period_us").write_text("100000")
    proc = proc_fixture(tmp_path, "4:memory:/container\n3:cpu,cpuacct:/\n",
        f"1 0 0:1 /container {memory} rw - cgroup cgroup rw,memory\n"
        f"2 0 0:2 /container {cpu} rw - cgroup cgroup rw,cpu,cpuacct\n")
    result = cgroup_limits(proc)
    assert result["memory_bytes"] == [2000000000]
    assert result["cpu_quotas"] == [1.0]
