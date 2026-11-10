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


def test_explicit_swap_changes_only_scheduler_not_tool_ram(tmp_path):
    ram = resolve(tmp_path / 'ram', cores=12, mem_mb=28672)
    swapped = resolve(tmp_path / 'swap', {'bactologix_compute': {'swap_mb': 98304}},
                      cores=12, mem_mb=126976)
    for name in ('total_mem_mb', 'assembly_mem_mb', 'spades_memory_gb',
                 'polishing_mem_mb', 'busco_mem_mb', 'prokka_mem_mb', 'kraken_mem_mb'):
        assert swapped.config[name] == ram.config[name]
    assert swapped.physical_memory_budget == ram.memory_budget == 25600
    assert swapped.memory_budget == 123904  # Keep native 80% physical RAM policy.
    assert swapped.record['memory_policy_at_resolution'] == {
        'requested_scheduler_mb': 126976, 'requested_ram_mb': 28672,
        'additional_swap_mb': 98304, 'effective_physical_ram_mb': 25600,
        'effective_scheduler_mb': 123904,
        'capacity_check': 'Native physical-RAM policy applied; usable swap must be checked by the launcher',
    }
    assert 'physical RAM budget=25600' in swapped.summary()
    assert 'additional swap allowance=98304' in swapped.summary()
    assert 'effective scheduler memory budget=123904' in swapped.summary()


def test_swap_does_not_suppress_explicit_physical_limits_or_spades_checks(tmp_path):
    values = {'bactologix_compute': {'swap_mb': 98304}, 'total_mem_mb': 16000,
              'assembly_mem_mb': 10000, 'spades_memory_gb': 7, 'polishing_mem_mb': 2000}
    resolved = resolve(tmp_path, values, mem_mb=126976)
    assert resolved.physical_memory_budget == 16000
    assert resolved.memory_budget == 114304
    assert resolved.config['assembly_mem_mb'] == 10000
    assert resolved.config['spades_memory_gb'] == 7
    assert resolved.config['polishing_mem_mb'] == 2000
    with pytest.raises(ValueError, match='detected host/container'):
        resolve(tmp_path, {**values, 'total_mem_mb': 40000}, mem_mb=126976)
    with pytest.raises(ValueError, match='fit within'):
        resolve(tmp_path, {**values, 'assembly_mem_mb': 1024}, mem_mb=126976)


@pytest.mark.parametrize('policy', [None, [], 'swap', {'swap_mb': -1}, {'swap_mb': True},
                                    {'swap_mb': '98304'}, {'swap_mb': 1.5}, {'swap_mb': None}])
def test_malformed_swap_policy_is_rejected(tmp_path, policy):
    with pytest.raises(ValueError, match='bactologix_compute'):
        resolve(tmp_path, {'bactologix_compute': policy}, mem_mb=126976)


@pytest.mark.parametrize('budget', [None, 98303, 98304])
def test_swap_requires_an_explicit_positive_ram_portion(tmp_path, budget):
    with pytest.raises(ValueError, match='positive RAM budget'):
        resolve(tmp_path, {'bactologix_compute': {'swap_mb': 98304}}, mem_mb=budget)


def test_zero_swap_preserves_ram_only_resolution(tmp_path):
    ram = resolve(tmp_path / 'ram', mem_mb=12000)
    zero = resolve(tmp_path / 'zero', {'bactologix_compute': {'swap_mb': 0}}, mem_mb=12000)
    assert ram.memory_budget == zero.memory_budget == zero.physical_memory_budget == 12000
    assert 'scheduler_memory_policy' not in zero.identity


def test_swap_resume_preserves_record_and_rechecks_ram_ceiling(tmp_path):
    values = {'bactologix_compute': {'swap_mb': 98304}}
    first = resolve(tmp_path, values, mem_mb=126976)
    first.persist()
    original = first.path.read_bytes()
    resumed = resolve(tmp_path, values, mem_mb=110304, capacity=capacity(memory_mb=16000))
    assert resumed.reused
    assert resumed.config == first.config
    assert resumed.physical_memory_budget == 12000
    assert resumed.memory_budget == 110304
    assert resumed.record['memory_policy_at_resolution'] == first.record['memory_policy_at_resolution']
    resumed.persist()
    assert first.path.read_bytes() == original
    assert resumed.identity['scheduler_memory_policy'] == 'physical_ram_plus_explicit_swap_v1'


def test_swap_worker_reuses_controller_receipt_instead_of_rule_resource(tmp_path):
    values = {'bactologix_compute': {'swap_mb': 98304}}
    with pytest.raises(ValueError, match='controller.*frozen settings'):
        resolve(tmp_path, values, mem_mb=1024, controller=False)
    first = resolve(tmp_path, values, mem_mb=126976)
    first.persist()
    frozen = first.path.read_bytes()
    worker = resolve(tmp_path, values, mem_mb=1024, controller=False)
    assert worker.reused and worker.ram_scheduler_mb == 28672
    assert worker.memory_budget == first.memory_budget
    assert worker.config == first.config
    assert first.path.read_bytes() == frozen


def test_worker_rejects_a_missing_or_inconsistent_swap_receipt(tmp_path):
    values = {'bactologix_compute': {'swap_mb': 98304}}
    first = resolve(tmp_path, values, mem_mb=126976)
    first.persist()
    record = json.loads(first.path.read_text())
    record.pop('record_sha256')
    record['memory_policy_at_resolution']['requested_ram_mb'] += 1
    first.path.write_text(json.dumps({**record, 'record_sha256': digest(record)}))
    with pytest.raises(ValueError, match='controller swap scheduling receipt'):
        resolve(tmp_path, values, mem_mb=1024, controller=False)


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
