"""Native BUSCO retention checks using synthetic I/O, never genome analysis.

HYBAS_BUSCO_SOURCE enables source-backed checks against the pinned native code.
"""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/lib"))
import busco_retention as adapter
from native_output_retention import RetentionError, RetentionSession


@pytest.fixture
def native_source():
    value = os.environ.get("HYBAS_BUSCO_SOURCE")
    if value:
        source = Path(value)
        for relative in ("src/busco", "busco"):
            if (source / relative / "_version.py").is_file():
                source = source / relative
                break
    else:
        try:
            source = adapter.source_directory()
        except adapter.UnsupportedBuscoError:
            pytest.skip("Set HYBAS_BUSCO_SOURCE or install BUSCO 6.0.0")
    adapter.verify_source(source)  # Explicit invalid source must fail, never skip.
    return source


@pytest.fixture
def prepared(tmp_path, native_source):
    root = tmp_path / "output"
    root.mkdir()
    session = RetentionSession.create(root, tool="busco")
    entry, env, path, receipt = adapter.prepare(native_source, session, ["-m", "genome"])
    yield SimpleNamespace(root=root, session=session, entry=entry, env=env,
                          receipt_path=path, receipt=receipt)
    session.close()


def run_source_script(prepared, code, *, name="check.py"):
    script = prepared.root.parent / name
    script.write_text(code)
    return subprocess.run([sys.executable, str(script), str(prepared.root)],
                          cwd=prepared.root.parent, env=prepared.env,
                          capture_output=True, text=True, timeout=30)


def require_busco_dependencies():
    for module in ("Bio", "numpy", "pandas", "requests"):
        if importlib.util.find_spec(module) is None:
            if os.environ.get("HYBAS_BUSCO_SOURCE"):
                pytest.fail(f"Explicit native BUSCO checks require dependency {module}")
            pytest.skip(f"BUSCO dependency {module} is needed")


def events(session):
    return [(path, json.loads((path / "manifest.json").read_text()))
            for path in (Path(session.path) / "events").iterdir()]


def test_native_call_arguments_and_result_object_are_unchanged():
    calls, captured = [], []
    result = subprocess.CompletedProcess(["native", "arg"], 0, b"out", b"err")

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        return result

    session = SimpleNamespace(capture_process=lambda command, value: captured.append((command, value)))
    proxy = adapter.CapturedSubprocess(SimpleNamespace(run=run, PIPE=subprocess.PIPE), session)
    command = ["native", "--threshold", "unchanged"]
    kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, cwd="original")
    assert proxy.run(command, **kwargs) is result
    assert calls == [((command,), kwargs)]
    assert captured == [(command, result)]
    assert proxy.PIPE == subprocess.PIPE


@pytest.mark.parametrize("error", [
    subprocess.CalledProcessError(7, ["native"], output=b"failed out", stderr=b"failed err"),
    subprocess.TimeoutExpired(["native"], 0.1, output=b"partial out", stderr=b"partial err"),
])
def test_native_errors_and_partial_streams_keep_the_same_exception(error):
    captured = []

    def run(*args, **kwargs):
        raise error

    proxy = adapter.CapturedSubprocess(SimpleNamespace(run=run),
        SimpleNamespace(capture_process=lambda command, value: captured.append((command, value))))
    with pytest.raises(type(error)) as raised:
        proxy.run(["native"], check=True)
    assert raised.value is error
    assert captured == [(["native"], error)]


def test_wrong_version_fails_before_import_or_launch(tmp_path):
    (tmp_path / "_version.py").write_text('__version__ = "999.0"\n')
    with pytest.raises(adapter.UnsupportedBuscoError, match="Expected BUSCO"):
        adapter.verify_source(tmp_path)


def test_source_guard_and_private_copy_leave_installation_unchanged(prepared, native_source):
    assert prepared.receipt["busco_version"] == "6.0.0"
    package = Path(prepared.session.path) / "driver/busco"
    for relative, expected in adapter.REVIEWED_SOURCES.items():
        original = (native_source / relative).read_bytes()
        assert adapter.sha256(native_source / relative) == expected
        assert (package / relative).read_bytes() == original + adapter.FOOTER.encode()
        assert prepared.receipt["source_files"][relative]["source_sha256"] == expected
    assert prepared.receipt["adapter_files"]["busco_retention.py"] == adapter.sha256(adapter.__file__)
    assert (package / "ConfigManager.py").read_bytes() == (native_source / "ConfigManager.py").read_bytes()


@pytest.mark.parametrize("module", ["busco_tools/Toolset.py", "run_BUSCO.py"])
def test_changed_native_cleanup_source_is_rejected(tmp_path, native_source, module):
    source = tmp_path / "busco"
    shutil.copytree(native_source, source)
    target = source / module
    target.write_bytes(target.read_bytes() + b"\n# unreviewed change\n")
    with pytest.raises(adapter.UnsupportedBuscoError, match="Unreviewed BUSCO source"):
        adapter.verify_source(source)


def test_additional_native_module_is_rejected(tmp_path, native_source):
    source = tmp_path / "busco"
    shutil.copytree(native_source, source)
    (source / "extra_cleanup.py").write_text("# unreviewed module\n")
    with pytest.raises(adapter.UnsupportedBuscoError, match="additional unreviewed"):
        adapter.verify_source(source)


def test_original_native_cleanup_paths_disappear_but_bytes_survive(prepared):
    require_busco_dependencies()
    result = run_source_script(prepared, r'''
from pathlib import Path
from types import SimpleNamespace
import sys
from busco.analysis.BuscoAnalysis import BuscoAnalysis
from busco.BuscoConfig import BuscoConfigMain
from busco.busco_tools.base import BaseRunner, GenePredictor
from busco.busco_tools.augustus import AugustusRunner, GFF2GBRunner
from busco.busco_tools.metaeuk import MetaeukRunner
from busco.busco_tools.sepp import SEPPRunner
root = Path(sys.argv[1])
cases = [
    ('force', lambda p: BuscoConfigMain._force_remove_existing_output_dir(str(p))),
    ('recreate', lambda p: BaseRunner.create_dirs(str(p), overwrite=True)),
    ('augustus', lambda p: AugustusRunner.remove_individual_err_logs(SimpleNamespace(tmp_dir=str(p)))),
    ('gff2gb', lambda p: GFF2GBRunner.remove_tmp_dir(SimpleNamespace(tmp_dir=str(p)))),
    ('metaeuk', lambda p: MetaeukRunner.remove_tmp_files(SimpleNamespace(_tmp_folder=str(p)))),
    ('sepp', lambda p: SEPPRunner.cleanup(SimpleNamespace(_tmp_folder=str(p)))),
]
for name, cleanup in cases:
    folder = root / name
    folder.mkdir()
    (folder / 'product.txt').write_text('native-' + name)
    cleanup(folder)
    assert not (folder / 'product.txt').exists()
    assert folder.is_dir() == (name == 'recreate')
main = root / 'run'
temp = main / 'tmp_dir'
temp.mkdir(parents=True)
(temp / 'product.txt').write_text('native-common-cleanup')
lineages = root / 'database/lineages'
lineages.mkdir(parents=True)
(lineages / 'old-download.tar').write_text('database archive stays untouched')
BuscoAnalysis.cleanup(SimpleNamespace(main_out=str(main), _lineages_download_path=str(lineages)))
assert not temp.exists()
assert (lineages / 'old-download.tar').read_text() == 'database archive stays untouched'
reference = root / 'reference.faa'
reference.write_text('reference payload')
temp = root / 'tmp'
temp.mkdir()
(temp / reference.name).write_text('previous-decompressed-reference')
returned = GenePredictor.decompress_refseq_file(SimpleNamespace(tmp_dir=str(temp)), str(reference))
assert Path(returned).is_symlink()
assert Path(returned).read_text() == 'reference payload'
''')
    assert result.returncode == 0, result.stdout + result.stderr
    prepared.session.validate()
    saved = {path.read_text() for path in Path(prepared.session.path).glob("events/*/files/**/*") if path.is_file() and not path.is_symlink()}
    assert {"native-" + name for name in ("force", "recreate", "augustus", "gff2gb", "metaeuk", "sepp", "common-cleanup")} <= saved
    assert "previous-decompressed-reference" in saved


def test_spawned_original_jobs_keep_unique_success_and_failed_streams(prepared):
    (prepared.root / "shared-out.log").write_text("older stdout baseline")
    (prepared.root / "shared-err.log").write_text("older stderr baseline")
    result = run_source_script(prepared, r'''
from pathlib import Path
from multiprocessing import Pool, Value
import sys
from busco.busco_tools.Toolset import Job, Tool
from busco.Exceptions import BuscoError

def worker(task):
    number, root = task
    command = "import sys; print('stdout-%d'); print('stderr-%d', file=sys.stderr); sys.exit(%d)" % (number, number, 7 if number == 3 else 0)
    job = Job('synthetic', sys.executable, str(Path(root) / 'shared-out.log'), str(Path(root) / 'shared-err.log'), None, root)
    job.add_parameter('-c')
    job.add_parameter(command)
    try:
        job.run()
        return 'success'
    except BuscoError:
        return 'failed'

if __name__ == '__main__':
    with Pool(2, initializer=Tool.init_globals, initargs=(Value('i', 0),)) as pool:
        statuses = pool.map(worker, [(number, sys.argv[1]) for number in range(4)])
    assert statuses == ['success', 'success', 'success', 'failed']
''', name="spawn-check.py")
    assert result.returncode == 0, result.stdout + result.stderr
    prepared.session.validate()
    captures = [(path, data) for path, data in events(prepared.session) if data["operation"] == "native-process-streams"]
    assert len(captures) == 4
    assert 1 <= len({data["pid"] for _, data in captures}) <= 2
    assert sorted(data["returncode"] for _, data in captures) == [0, 0, 0, 7]
    assert {path.joinpath("stdout.bin").read_bytes() for path, _ in captures} == {f"stdout-{i}\n".encode() for i in range(4)}
    assert {path.joinpath("stderr.bin").read_bytes() for path, _ in captures} == {f"stderr-{i}\n".encode() for i in range(4)}
    assert all(data["command"][0] == sys.executable for _, data in captures)
    baseline_copies = [(path / data["archive"]).read_text()
                      for path, data in events(prepared.session) if data.get("kind") == "file"]
    assert sorted(baseline_copies) == ["older stderr baseline", "older stdout baseline"]


def test_original_log_move_preserves_previous_destination(prepared):
    require_busco_dependencies()
    result = run_source_script(prepared, r'''
from pathlib import Path
from configparser import ConfigParser
import sys
from busco.BuscoRunner import AnalysisRunner
from busco.BuscoLogger import BuscoLogger
root = Path(sys.argv[1])
folder = root / 'run/logs'
folder.mkdir(parents=True)
(folder / 'busco.log').write_text('previous-busco-log')
source = Path('busco_' + BuscoLogger.pid + '.log')
source.write_text('current-busco-log')
config = ConfigParser()
config['busco_run'] = {'out_path': str(root), 'out': 'run'}
AnalysisRunner.move_log_file(config)
assert not source.exists()
assert (folder / 'busco.log').read_text() == 'current-busco-log'
''')
    assert result.returncode == 0, result.stdout + result.stderr
    prepared.session.validate()
    assert any((path / data["archive"]).read_text() == "previous-busco-log" for path, data in events(prepared.session) if data.get("kind") == "file")


def test_native_cleanup_archive_failure_keeps_product_and_latches_failure(prepared):
    require_busco_dependencies()
    result = run_source_script(prepared, r'''
from pathlib import Path
import sys
import native_output_retention as retention
from busco.busco_tools.base import BaseRunner
folder = Path(sys.argv[1]) / 'native-output'
folder.mkdir()
(folder / 'product').write_text('must remain')
def unavailable_archive(*args, **kwargs):
    raise OSError('synthetic disk full')
retention._write_json = unavailable_archive
try:
    BaseRunner.create_dirs(str(folder), overwrite=True)
except retention.RetentionError:
    pass
else:
    raise AssertionError('native OSError handler swallowed preservation failure')
assert (folder / 'product').read_text() == 'must remain'
''')
    assert result.returncode == 0, result.stdout + result.stderr
    with pytest.raises(RetentionError, match="preservation failed"):
        prepared.session.validate()


def test_native_timeout_preserves_partial_output(prepared):
    result = run_source_script(prepared, r'''
from pathlib import Path
from multiprocessing import Value
import sys
from busco.busco_tools.Toolset import Job, Tool
root = Path(sys.argv[1])
Tool.init_globals(Value('i', 0))
job = Job('synthetic', sys.executable, str(root / 'timeout-out.log'), str(root / 'timeout-err.log'), 0.3, str(root))
job.add_parameter('-c')
job.add_parameter("import sys,time; print('partial-out',flush=True); print('partial-err',file=sys.stderr,flush=True); time.sleep(3)")
job.run()
''')
    assert result.returncode == 0, result.stdout + result.stderr
    prepared.session.validate()
    captures = [(path, data) for path, data in events(prepared.session) if data["operation"] == "native-process-streams"]
    assert len(captures) == 1
    path, receipt = captures[0]
    assert receipt["outcome"] == "TimeoutExpired"
    assert receipt["partial"] is True
    assert (path / "stdout.bin").read_bytes() == b"partial-out\n"
    assert (path / "stderr.bin").read_bytes() == b"partial-err\n"


def test_launcher_preserves_arguments_returncode_and_checks_failure_latch(tmp_path, native_source):
    root = tmp_path / "output"
    root.mkdir()
    arguments = ["-i", "relative-input.fa", "-m", "genome", "-o", "run", "--out_path", str(root), "--offline", "--force"]
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=9)

    assert adapter.launch(arguments, root, source=native_source, runner=runner) == 9
    assert calls[0][0][2:] == arguments
    assert "cwd" not in calls[0][1]
    assert json.loads(next(root.glob("native-retained/*/busco-driver.json")).read_text())["status"] == "failed"

    def caught_failure(command, env):
        session = RetentionSession.attach(env["HYBAS_NATIVE_RETENTION_SESSION"])
        try:
            session.remove(tmp_path / "external.txt")
        except RetentionError:
            pass  # Simulate upstream swallowing a preservation failure.
        finally:
            session.close()
        return SimpleNamespace(returncode=0)

    external = tmp_path / "external.txt"
    external.write_text("must survive")
    with pytest.raises(RetentionError, match="preservation failed"):
        adapter.launch(arguments, root, source=native_source, runner=caught_failure)
    assert external.read_text() == "must survive"


@pytest.mark.parametrize("name", ["..", ".", "native-retained", "native-retained/sub"])
def test_launcher_rejects_output_that_could_remove_control_or_escape_root(tmp_path, name):
    with pytest.raises(ValueError, match="separate directory"):
        adapter.launch(["-o", name, "--out_path", str(tmp_path)], tmp_path)


def test_launcher_rejects_symlinked_output_root_before_mutation(tmp_path):
    owned = tmp_path / "real"
    owned.mkdir()
    link = tmp_path / "link"
    link.symlink_to(owned, target_is_directory=True)
    with pytest.raises((ValueError, RetentionError)):
        adapter.launch(["-o", "run", "--out_path", str(link)], link)
    assert not list(owned.iterdir())
