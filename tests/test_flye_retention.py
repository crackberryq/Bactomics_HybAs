"""Retention-only contracts; native source tests use synthetic, stubbed analysis.

Set HYBAS_FLYE_SOURCE to an official Flye 2.9.6 checkout (or an installed
package's parent). No source is downloaded and no biological tool is executed.
"""
from contextlib import nullcontext
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/lib"))
import flye_retention as adapter
import native_output_retention as retention
from native_output_retention import RetentionError, RetentionSession


def files(root):
    return {path.relative_to(root).as_posix(): path.read_bytes()
            for path in root.rglob("*")
            if path.is_file() and "native-retained" not in path.relative_to(root).parts}


def assert_saved(session, contents):
    # Check payload bytes independently of the helper's receipt serialization.
    assert any(path.is_file() and path.read_bytes() == contents
               for path in session.path.rglob("*"))


@pytest.fixture
def native(monkeypatch):
    source = os.environ.get("HYBAS_FLYE_SOURCE")
    if not source:
        pytest.skip("Set HYBAS_FLYE_SOURCE for version/hash-bound upstream source tests")
    monkeypatch.syspath_prepend(source)
    return adapter.verified_modules()


def test_rejects_version_before_loading_native_workflow(monkeypatch):
    seen = []

    def load(name):
        seen.append(name)
        return SimpleNamespace(__version__="2.9.7")

    monkeypatch.setattr(adapter.importlib, "import_module", load)
    with pytest.raises(RetentionError, match="reviewed version"):
        adapter.verified_modules()
    assert seen == ["flye.__version__"]


def test_rejects_changed_native_source_before_entry(tmp_path, monkeypatch):
    source = tmp_path / "main.py"
    source.write_text("# unreviewed native source\n")
    root = tmp_path / "output"
    root.mkdir()
    monkeypatch.setattr(adapter.importlib, "import_module", lambda name:
                        SimpleNamespace(__version__="2.9.6", __file__=str(source)))
    with pytest.raises(RetentionError, match="source guard failed"):
        adapter.run(root, ["--out-dir", str(root)])
    assert not (root / "native-retained").exists()


@pytest.mark.parametrize("reason", ["occupied", "mismatch", "symlink"])
def test_rejects_nonfresh_or_mismatched_root_before_native_import(tmp_path, monkeypatch, reason):
    root = tmp_path / "output"
    root.mkdir()
    output = root
    if reason == "occupied":
        (root / "old.fasta").write_bytes(b"old result")
    elif reason == "mismatch":
        output = tmp_path / "different"
    else:
        link = tmp_path / "link"
        link.symlink_to(root, target_is_directory=True)
        root = output = link
    monkeypatch.setattr(adapter, "verified_modules", lambda:
                        pytest.fail("Native modules must not load for an invalid output root"))
    with pytest.raises(RetentionError):
        adapter.run(root, ["--out-dir", str(output)])
    assert not (root / "native-retained").exists()


def test_native_arguments_and_module_bindings_are_restored(tmp_path, monkeypatch):
    native_args = ["--pacbio-hifi", "reads with spaces.fq", "--out-dir", str(tmp_path),
                   "--threads", "12", "--iterations", "1"]
    old_argv = sys.argv
    modules = {name: SimpleNamespace(os=os) for name in adapter.MODULE_SHA256}

    def entry():
        assert sys.argv == ["flye", *native_args]
        assert modules["flye.main"].os is not os
        assert modules["flye.polishing.alignment"].os is os
        return 7

    modules["flye.main"].main = entry
    monkeypatch.setattr(adapter, "verified_modules", lambda: modules)
    assert adapter.run(tmp_path, native_args) == 7
    assert sys.argv is old_argv
    assert all(module.os is os for module in modules.values())
    assert all(not hasattr(module, "open") for module in modules.values())
    receipt = next((tmp_path / "native-retained").rglob("session.json"))
    metadata = json.loads(receipt.read_text())["metadata"]
    assert metadata["native_argv"] == ["flye", *native_args]
    assert metadata["native_module_sha256"] == adapter.MODULE_SHA256
    assert metadata["native_version"] == "2.9.6"
    assert metadata["adapter_sha256"] == hashlib.sha256(Path(adapter.__file__).read_bytes()).hexdigest()


def test_upstream_catch_cannot_hide_preservation_failure(tmp_path, monkeypatch):
    modules = {name: SimpleNamespace(os=os) for name in adapter.MODULE_SHA256}
    outside = tmp_path.parent / (tmp_path.name + "-external")
    outside.write_bytes(b"untouched external data")
    old_argv = sys.argv

    def entry():
        try:
            modules["flye.main"].os.remove(outside)
        except Exception:
            pass  # Simulate a broad upstream cleanup catch.
        return 0

    modules["flye.main"].main = entry
    monkeypatch.setattr(adapter, "verified_modules", lambda: modules)
    with pytest.raises(RetentionError):
        adapter.run(tmp_path, ["--out-dir", str(tmp_path)])
    assert outside.read_bytes() == b"untouched external data"
    assert sys.argv is old_argv
    assert all(module.os is os for module in modules.values())


def test_native_consensus_cleanup_preserves_bytes_and_original_outputs(native, tmp_path, monkeypatch):
    main, aln = native["flye.main"], native["flye.polishing.alignment"]
    draft = tmp_path / "input.fasta"
    draft.write_text(">contig_1\nACGT\n")
    bam = b"synthetic consensus mapping\x00\xff"

    def mapping(reference, reads, threads, platform, read_type, output):
        Path(output).write_bytes(bam)
        Path(output + ".csi").write_bytes(b"synthetic index")

    monkeypatch.setattr(aln, "make_alignment", mapping)
    monkeypatch.setattr(main.cons, "get_consensus", lambda *args: {"contig_1": "ACGT"})
    args = SimpleNamespace(reads=["synthetic.fastq"], threads=1, platform="pacbio", read_type="hifi")
    roots = [tmp_path / "control", tmp_path / "retained"]
    for root in roots:
        root.mkdir()
    session = RetentionSession.create(roots[1], tool="flye")
    for root, context in zip(roots, [nullcontext(), adapter.preserve_native_outputs(native, session)]):
        with context:
            main.JobConsensus(args, str(root), str(draft)).run()
        assert not (root / "10-consensus/minimap.bam").exists()
    assert files(roots[0]) == files(roots[1])
    assert_saved(session, bam)


def test_native_polishing_cleanup_preserves_every_reviewed_product(native, tmp_path, monkeypatch):
    main, pol = native["flye.main"], native["flye.polishing.polish"]
    draft, edges, graph = [tmp_path / name for name in ("input.fasta", "edges.fasta", "graph.gfa")]
    draft.write_text(">contig_1\nACGT\n>contig_2\nTT\n")
    edges.write_text(">edge_1\nACGT\n")
    graph.write_text("H\tVN:Z:1.0\nS\tedge_1\tACGT\tdp:i:4\n")
    generated = {}

    def mapping(reference, reads, threads, platform, read_type, output):
        content = b"synthetic mapping: " + Path(output).name.encode()
        Path(output).write_bytes(content)
        Path(output + ".csi").write_bytes(b"index: " + Path(output).name.encode())
        generated[Path(output).name] = content

    def bubbles(alignment, info, assembly, platform, threads, output):
        content = b"synthetic bubbles: " + Path(output).name.encode()
        Path(output).write_bytes(content)
        generated[Path(output).name] = content
        return {"contig_1": 100, "contig_2": 1}, 0.0

    def polish_bin(bubbles, matrix, hopo, output, threads, progress, use_hopo):
        content = b"synthetic consensus: " + Path(output).name.encode()
        Path(output).write_bytes(content)
        generated[Path(output).name] = content

    monkeypatch.setattr(pol, "make_alignment", mapping)
    monkeypatch.setattr(pol, "make_bubbles", bubbles)
    monkeypatch.setattr(pol, "_run_polish_bin", polish_bin)
    monkeypatch.setattr(pol, "_compose_sequence", lambda output:
                        ({"contig_1": "ACGT", "contig_2": "TT"}, {"contig_1": 4, "contig_2": 2},
                         {"contig_1": [(4, 100)], "contig_2": [(2, 1)]}))
    monkeypatch.setattr(pol, "SynchronizedSamReader", lambda *args:
                        SimpleNamespace(get_all_alignments=lambda: []))
    monkeypatch.setattr(pol, "multiprocessing", SimpleNamespace(Manager=lambda: object()))
    # Freeze gzip's header timestamp so the complete final files can be compared.
    monkeypatch.setattr(gzip, "time", SimpleNamespace(time=lambda: 1234567890))
    args = SimpleNamespace(reads=["synthetic.fastq"], threads=1, platform="pacbio",
                           read_type="hifi", num_iters=2, meta=False)
    roots = [tmp_path / "control", tmp_path / "retained"]
    for root in roots:
        root.mkdir()
    session = RetentionSession.create(roots[1], tool="flye")
    for root, context in zip(roots, [nullcontext(), adapter.preserve_native_outputs(native, session)]):
        with context:
            main.JobPolishing(args, str(root), "unused.log", str(draft), str(edges), str(graph)).run()
    assert files(roots[0]) == files(roots[1])
    assert len(generated) == 7  # Two iteration BAM/bubbles/consensus sets + graph mapping.
    for name, content in generated.items():
        assert not (roots[1] / "40-polishing" / name).exists()
        assert_saved(session, content)
    assert not (roots[1] / "40-polishing/polished_2.fasta").exists()
    assert_saved(session, draft.read_bytes())  # Last unfiltered contigs include the removed contig.
    assert b"contig_2" not in (roots[1] / "40-polishing/filtered_contigs.fasta").read_bytes()


def test_native_mapping_stderr_revisions_and_command_equivalence(native, tmp_path, monkeypatch):
    aln = native["flye.polishing.alignment"]
    reference = tmp_path / "reference.fasta"
    reference.write_text(">contig_1\nACGT\n")
    commands = []
    count = 0

    def check_call(command, **kwargs):
        nonlocal count
        commands.append(command)
        if "stderr" in kwargs:
            count += 1
            kwargs["stderr"].write(f"native mapping log {count}\n")
            kwargs["stderr"].flush()

    monkeypatch.setattr(aln, "subprocess", SimpleNamespace(check_call=check_call,
                        CalledProcessError=subprocess.CalledProcessError))
    monkeypatch.setattr(aln, "datetime", SimpleNamespace(datetime=SimpleNamespace(
                        now=lambda: SimpleNamespace(strftime=lambda pattern: "fixed-sort-prefix"))))
    roots = [tmp_path / "control", tmp_path / "retained"]
    for root in roots:
        root.mkdir()
    session = RetentionSession.create(roots[1], tool="flye")
    for root, context in zip(roots, [nullcontext(), adapter.preserve_native_outputs(native, session)]):
        count = 0
        with context:
            for index in (1, 2):
                aln._run_minimap(str(reference), ["synthetic.fastq"], 2, "pacbio-hifi",
                                 str(root / f"minimap_{index}.bam"))
    assert files(roots[0]) == files(roots[1])
    assert_saved(session, b"native mapping log 1\n")
    normalized = [str(command).replace(str(roots[0]), "OUT").replace(str(roots[1]), "OUT")
                  for command in commands]
    assert normalized[:4] == normalized[4:]


def test_native_state_revisions_and_repeated_cleanup_are_distinct(native, tmp_path, monkeypatch):
    main = native["flye.main"]
    monkeypatch.setattr(main.Job, "run_params", {"stage_name": "", "synthetic": True})
    session = RetentionSession.create(tmp_path, tool="flye")
    prior_states = []
    with adapter.preserve_native_outputs(native, session):
        for stage in ("configure", "consensus", "polishing"):
            job = main.Job()
            job.name = stage
            job.save(str(tmp_path / "params.json"))
            prior_states.append((tmp_path / "params.json").read_bytes())
        for content in (b"first mapping", b"second mapping"):
            path = tmp_path / "minimap.bam"
            path.write_bytes(content)
            main.os.remove(path)
            assert not path.exists()
    for content in prior_states[:-1] + [b"first mapping", b"second mapping"]:
        assert_saved(session, content)
    assert json.loads((tmp_path / "params.json").read_text())["stage_name"] == "polishing"


def test_native_failure_restores_bindings_without_deleting_partial_results(native, tmp_path, monkeypatch):
    main = native["flye.main"]
    original_argv = sys.argv
    original_os = main.os

    def fail():
        (tmp_path / "partial.fasta").write_bytes(b"partial native output")
        raise RuntimeError("synthetic native failure")

    monkeypatch.setattr(main, "main", fail)
    with pytest.raises(RuntimeError, match="synthetic native failure"):
        adapter.run(tmp_path, ["--out-dir", str(tmp_path)])
    assert (tmp_path / "partial.fasta").read_bytes() == b"partial native output"
    assert main.os is original_os
    assert sys.argv is original_argv


@pytest.mark.parametrize("module_name", list(adapter.MODULE_SHA256))
def test_each_native_source_guard_is_enforced(native, tmp_path, monkeypatch, module_name):
    changed = tmp_path / "changed.py"
    changed.write_text("# installation changed after review\n")
    monkeypatch.setattr(native[module_name], "__file__", str(changed))
    root = tmp_path / "output"
    root.mkdir()
    with pytest.raises(RetentionError, match=module_name):
        adapter.run(root, ["--out-dir", str(root)])
    assert not list(root.iterdir())


def test_archive_io_failure_keeps_original_and_fails_native_invocation(native, tmp_path, monkeypatch):
    main = native["flye.main"]
    output = tmp_path / "minimap.bam"
    original_rename = retention.os.rename

    def fail_archive(source, target, *args, **kwargs):
        if Path(source) == output:
            raise OSError("synthetic full disk")
        return original_rename(source, target, *args, **kwargs)

    def entry():
        output.write_bytes(b"native alignment that must survive")
        try:
            main.os.remove(output)
        except Exception:
            pass
        return 0

    monkeypatch.setattr(main, "main", entry)
    monkeypatch.setattr(retention.os, "rename", fail_archive)
    with pytest.raises(RetentionError):
        adapter.run(tmp_path, ["--out-dir", str(tmp_path)])
    assert output.read_bytes() == b"native alignment that must survive"
