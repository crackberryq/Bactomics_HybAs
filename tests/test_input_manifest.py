"""Metadata-only v10 input tests; no biological sequences or tool invocations."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
from input_manifest import ManifestError, MODES, load_manifest, validate_manifest


def library(name="pe1", platform="illumina", **fields):
    group = {"read_group_id": "rg1", "run_id": "run1"}
    value = {"library_id": name, "platform": platform, "read_groups": [group]}
    if platform == "illumina":
        value.update(library_type="paired_end", orientation="fr")
        group.update(r1=f"reads/{name}_R1.fastq.gz", r2=f"reads/{name}_R2.fastq.gz")
    else:
        group["reads"] = f"reads/{name}.fastq.gz"
        if platform == "pacbio":
            value["chemistry"] = "hifi"
    value.update(fields)
    return value


def manifest(libraries=None, mode="auto"):
    return {"schema_version": 1, "samples": [{"sample_id": "example", "mode": mode,
            "libraries": libraries if libraries is not None else [library()]}]}


@pytest.mark.parametrize("mode,libraries,chemistry", [
    ("illumina_only", [library()], None),
    ("dual_illumina_runs", [library(), library("pe2")], None),
    ("dual_illumina_complementary", [library(insert_size_bp=350),
        library("mp", library_type="mate_pair", orientation="rf", insert_size_bp=3000)], None),
    ("dual_illumina_complementary", [library(insert_size_bp=350), library("pe2", insert_size_bp=700)], None),
    ("ont_only", [library("ont", "ont")], None),
    ("illumina_ont", [library(), library("ont", "ont")], None),
    ("pacbio_only", [library("pb", "pacbio")], "hifi"),
    ("pacbio_only", [library("pb", "pacbio", chemistry="clr")], "clr"),
    ("illumina_pacbio", [library(), library("pb", "pacbio")], "hifi"),
    ("illumina_pacbio", [library(), library("pb", "pacbio", chemistry="clr")], "clr"),
])
def test_supported_metadata_routes(tmp_path, mode, libraries, chemistry):
    source = manifest(libraries, mode)
    before = deepcopy(source)
    result = validate_manifest(source, base_dir=tmp_path)
    selection = result["selections"][0]
    assert selection["selected_mode"] == mode
    assert selection["pacbio_chemistry"] == chemistry
    assert selection["library_ids"] == [x["library_id"] for x in libraries]
    assert selection["execution_ready"] is result["execution_ready"] is False
    assert selection["qualification"] == "manifest_validation_only"
    assert source == before
    assert list(tmp_path.iterdir()) == []  # No reads or workflow outputs needed/created.


def test_lanes_and_runs_are_not_confused(tmp_path):
    pe = library()
    second = dict(pe["read_groups"][0], read_group_id="rg2", lane_id="L002",
                  r1="reads/lane2_R1.fq", r2="reads/lane2_R2.fq")
    pe["read_groups"].append(second)
    assert validate_manifest(manifest([pe]), base_dir=tmp_path)["selections"][0]["selected_mode"] == "illumina_only"
    second["run_id"] = "run2"
    result = validate_manifest(manifest([pe]), base_dir=tmp_path)
    assert result["selections"][0]["selected_mode"] == "dual_illumina_runs"
    assert len(result["manifest"]["samples"][0]["libraries"]) == 1
    assert len(result["manifest"]["samples"][0]["libraries"][0]["read_groups"]) == 2


@pytest.mark.parametrize("libraries,expected", [
    ([library()], "illumina_only"),
    ([library("ont", "ont")], "ont_only"),
    ([library(), library("ont", "ont")], "illumina_ont"),
    ([library("pb", "pacbio")], "pacbio_only"),
    ([library(), library("pb", "pacbio")], "illumina_pacbio"),
])
def test_auto_uses_declarations(tmp_path, libraries, expected):
    result = validate_manifest(manifest(libraries), base_dir=tmp_path)
    assert result["selections"][0]["selected_mode"] == expected
    assert result["selections"][0]["selection_source"] == "declared_metadata"


@pytest.mark.parametrize("libraries,mode,message", [
    ([library(), library("pe2")], "auto", "explicit mode"),
    ([library()], "pacbio_only", "does not match"),
    ([library()], "dual_illumina_runs", "does not match"),
    ([library("ont", "ont"), library("pb", "pacbio")], "auto", "unsupported"),
    ([library(), library("ont", "ont"), library("pb", "pacbio")], "auto", "unsupported"),
    ([library("pb", "pacbio"), library("clr", "pacbio", chemistry="clr")], "auto", "mixed PacBio"),
    ([library(library_type="mate_pair")], "auto", "unsupported"),
    ([library(), library("mp", library_type="mate_pair"), library("ont", "ont")], "auto", "unsupported"),
    ([library(), library("pe2")], "dual_illumina_complementary", "insert_size_bp"),
    ([library(insert_size_bp=350), library("pe2", insert_size_bp=350)], "dual_illumina_complementary", "distinct"),
])
def test_incompatible_or_ambiguous_routes_fail(tmp_path, libraries, mode, message):
    with pytest.raises(ManifestError, match=message):
        validate_manifest(manifest(libraries, mode), base_dir=tmp_path)


@pytest.mark.parametrize("value", [True, False, 0, 2, "1", 1.0, None])
def test_schema_version_is_strict(tmp_path, value):
    data = manifest()
    data["schema_version"] = value
    with pytest.raises(ManifestError, match="schema_version"):
        validate_manifest(data, base_dir=tmp_path)


@pytest.mark.parametrize("field,value", [
    ("platform", "nanopore"), ("library_type", "single"), ("orientation", "unknown"),
    ("insert_size_bp", True), ("insert_size_bp", 0), ("insert_size_bp", -1),
    ("insert_size_bp", 350.0), ("chemistry", "hifi"), ("unexpected", 1),
])
def test_invalid_library_fields(tmp_path, field, value):
    lib = library()
    lib[field] = value
    with pytest.raises(ManifestError):
        validate_manifest(manifest([lib]), base_dir=tmp_path)


@pytest.mark.parametrize("value", [None, "auto", "HiFi", "", False])
def test_pacbio_chemistry_must_be_explicit(tmp_path, value):
    with pytest.raises(ManifestError, match="chemistry"):
        validate_manifest(manifest([library("pb", "pacbio", chemistry=value)]), base_dir=tmp_path)


@pytest.mark.parametrize("path", ["", " foo.fq", "file.fq\n", "reads/*.fq", "https://host/data.fq",
                                    "C:\\reads\\r1.fq", "C:/reads/r1.fq", "reads/ccs.bam", "reads/file.fa"])
def test_invalid_paths(tmp_path, path):
    lib = library()
    lib["read_groups"][0]["r1"] = path
    with pytest.raises(ManifestError):
        validate_manifest(manifest([lib]), base_dir=tmp_path)


def test_duplicate_files_across_samples_and_path_aliases(tmp_path):
    data = manifest()
    other = deepcopy(data["samples"][0])
    other["sample_id"] = "other"
    other["libraries"][0]["read_groups"][0]["r1"] = "reads/../reads/pe1_R1.fastq.gz"
    data["samples"].append(other)
    with pytest.raises(ManifestError, match="already declared"):
        validate_manifest(data, base_dir=tmp_path)


def test_symlink_alias_is_not_counted_twice(tmp_path):
    (tmp_path / "reads").mkdir()
    target = tmp_path / "reads/pe1_R1.fastq.gz"
    target.touch()
    (tmp_path / "alias.fq").symlink_to(target)
    lib = library()
    lib["read_groups"][0]["r2"] = "alias.fq"
    with pytest.raises(ManifestError, match="already declared"):
        validate_manifest(manifest([lib]), base_dir=tmp_path)


def test_unresolvable_path_is_manifest_error(tmp_path, monkeypatch):
    def fail_expanduser(_path):
        raise RuntimeError("Could not determine home directory")
    monkeypatch.setattr(Path, "expanduser", fail_expanduser)
    with pytest.raises(ManifestError, match="cannot resolve input path"):
        validate_manifest(manifest(), base_dir=tmp_path)


@pytest.mark.parametrize("field", ["r1", "r2", "run_id", "read_group_id"])
def test_missing_read_group_metadata(tmp_path, field):
    lib = library()
    del lib["read_groups"][0][field]
    with pytest.raises(ManifestError, match="missing fields"):
        validate_manifest(manifest([lib]), base_dir=tmp_path)


def test_long_read_library_rejects_paired_fields(tmp_path):
    with pytest.raises(ManifestError, match="only valid for Illumina"):
        validate_manifest(manifest([library("pb", "pacbio", orientation="fr")]), base_dir=tmp_path)


def test_pacbio_cannot_omit_chemistry(tmp_path):
    lib = library("pb", "pacbio")
    del lib["chemistry"]
    with pytest.raises(ManifestError, match="chemistry"):
        validate_manifest(manifest([lib]), base_dir=tmp_path)


@pytest.mark.parametrize("kind", ["sample", "library", "group", "mate"])
def test_duplicate_identifiers_and_mates(tmp_path, kind):
    data = manifest()
    sample = data["samples"][0]
    lib = sample["libraries"][0]
    if kind == "sample":
        data["samples"].append(deepcopy(sample))
    elif kind == "library":
        sample["libraries"].append(library("pe2", library_id="pe1"))
    elif kind == "group":
        lib["read_groups"].append(dict(lib["read_groups"][0]))
    else:
        lib["read_groups"][0]["r2"] = lib["read_groups"][0]["r1"]
    with pytest.raises(ManifestError, match="duplicate|already declared"):
        validate_manifest(data, base_dir=tmp_path)


@pytest.mark.parametrize("data", [None, [], {}, {"schema_version": 1, "samples": []},
    {"schema_version": 1, "samples": [None]}, {"schema_version": 1, "samples": [{"sample_id": "a", "libraries": []}]}])
def test_invalid_structures(tmp_path, data):
    with pytest.raises(ManifestError):
        validate_manifest(data, base_dir=tmp_path)


@pytest.mark.parametrize("suffix", [".json", ".yaml"])
def test_load_reproducibility_and_cli(tmp_path, suffix):
    path = tmp_path / ("inputs" + suffix)
    source = manifest()
    text = json.dumps(source) if suffix == ".json" else yaml.safe_dump(source)
    path.write_text(text)
    expected = load_manifest(path)
    assert expected["source"]["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert expected["manifest"]["samples"][0]["libraries"][0]["read_groups"][0]["r1"] == str(tmp_path / "reads/pe1_R1.fastq.gz")
    result = subprocess.run([sys.executable, str(ROOT / "workflow/scripts/validate_input_manifest.py"), str(path)],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == expected
    assert sorted(p.name for p in tmp_path.iterdir()) == [path.name]


@pytest.mark.parametrize("filename,text", [
    ("x.json", '{"schema_version": 1, "schema_version": 1, "samples": []}'),
    ("x.yaml", 'schema_version: 1\nschema_version: 1\nsamples: []\n'),
    ("x.yaml", 'schema_version: 1\nsamples: []\nfalse: bad\n'),
    ("x.json", '{bad'), ("x.yaml", 'samples: [}'),
])
def test_parse_errors_are_actionable(tmp_path, filename, text):
    path = tmp_path / filename
    path.write_text(text)
    result = subprocess.run([sys.executable, str(ROOT / "workflow/scripts/validate_input_manifest.py"), str(path)],
                            text=True, capture_output=True)
    assert result.returncode == 2
    assert "Input manifest error:" in result.stderr
    assert "Traceback" not in result.stderr
    assert not result.stdout


def test_shipped_example_covers_all_modes():
    result = load_manifest(ROOT / "examples/v10/inputs.example.yaml")
    assert {s["selected_mode"] for s in result["selections"]} == set(MODES)
    assert {s["pacbio_chemistry"] for s in result["selections"]} == {None, "hifi", "clr"}
