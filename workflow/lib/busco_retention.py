"""Run the reviewed BUSCO 6.0.0 driver with lossless output handling.

Only module-local filesystem and subprocess I/O are adapted. Commands, arguments,
result objects, scoring, and cleanup's absent-path semantics stay unchanged.
A private package copy applies the same adaptations to spawned workers. Native
programs' internal scratch files are outside this Python adapter's scope.
"""
from __future__ import annotations

import argparse
import ast
import builtins
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from native_output_retention import RetentionSession

BUSCO_VERSION = "6.0.0"
UPSTREAM_COMMIT = "6278721a1916f6da310e03ec9674099028c927a4"
# Identical in upstream 6.0.0 and Bioconda builds 0 and 3. Unadapted package
# files are copied unchanged and inventoried, including Conda's local config.
REVIEWED_SOURCES = {
    "Actions.py": "1ceca6a087058b8e4a70599d62393753526849000a0db3f2ffeb2fed5ef43590",
    "BuscoConfig.py": "1d5d3449374811c061b6e7c6245cf462b97d9ef6a2a4345121674a3012478b1f",
    "BuscoDownloadManager.py": "ac96110b499a0df688f52b893a50668fa64f4f5cee6c064466009b01dffe6afd",
    "BuscoPlacer.py": "f5cd98048c8f7dfdc8f986090945f00cd3572a567396c759526b35bf694f2d34",
    "BuscoRunner.py": "d707c318ae3ba021afdbd43ff955a70089972f6439475f1c5a7ef36e28405388",
    "analysis/BuscoAnalysis.py": "352382053d41876b7c64f3eb95919026dd529abf575bb6045f47945fc50718d0",
    "analysis/GenomeAnalysis.py": "b0b3418c080101f752a83760d058a79098ea5cf5f0738d4da61ed00f565c5b1d",
    "analysis/TranscriptomeAnalysis.py": "f4b81d6150192fa974b357fa27721e17d30b3bf3d3b0a4e7b55a1dc65126f996",
    "busco_tools/Toolset.py": "fa7599085ee4a0d2afa1896a138caef77ef0d83c8169011537cedd3dac430339",
    "busco_tools/augustus.py": "0bd346c2ffbedf85ff57972f5fa8eaca531e59dced6ac1a7a7986c5ecca4d5cd",
    "busco_tools/base.py": "54d6c073edae2d2d924fd38326ede776192bd13bf669f5078c90a56bc9d9e6ae",
    "busco_tools/blast.py": "6e5c8029236dd22f1ac8d719f105653b3033856d20297981bd233e5e1fbe0e9e",
    "busco_tools/hmmer.py": "8f3823379e9e464c0ca453e7fdcc976602259c7ffbc2df11173f2d0ae4d41531",
    "busco_tools/metaeuk.py": "3f2a79929378e872b7a926401b3f37895b1c6b884430f48f4ef3314c98894c2b",
    "busco_tools/miniprot.py": "2bdde54b489fa6b38ec4e419289617a8be437c5e033edf1ca8759766e8b44e8b",
    "busco_tools/prodigal.py": "d6fafe61af18049131e57dff2d27f40675dbbcb98ca71492c8348dd2cc2fbdfa",
    "busco_tools/sepp.py": "d7f23251407ab15dbcdd38ab061d20957e756d767ab261ee7a592caf8ad41e8a",
}
# Guard unadapted executable modules too: a new cleanup path must be reviewed.
SOURCE_HASHES = {name: (digest,) for name, digest in REVIEWED_SOURCES.items()}
SOURCE_HASHES.update({
    "AutoLineage.py": ("e8a4ba14b59ae159416364b860adb4bdb37a47fe3cf7ad6e6193f8dc7c925e56",),
    "BuscoLogger.py": ("c8b41f48b3b19e13f26a6f267bd38457f4bd7a5754f3ee7fbabdfc5cff7a9c29",),
    "BuscoPlot.py": ("975a5c62bc2af2bd018db4cb01621e7dbab7f962686d51e31ec67a990a329739",),
    "ConfigManager.py": ("121aff3194d1ba800ccaae1fae7a057837239e8be7751cab09ce68c688cbe418", "576f0138499eab75dc7c25a1d9715c0a80fc32762c9697b0b3c8ada0d5319636"),
    "Exceptions.py": ("bdbe602349baea502ae5f0f2594fc464d21b85548c79842ad778ef397654f7f9",),
    "__init__.py": ("4633517033a42295e5bb407cddb86f15b944a654ad6ae85c9c9bd44fba564b35",),
    "_version.py": ("be66f96d9b1449a56a12d8d526cddab1e0046d234e3a5be03560218bb94f4c13",),
    "analysis/Analysis.py": ("68652ee4b0bdab1763a533d79fcbba35458632acf0af71ad36f276859e8431dd",),
    "analysis/GeneSetAnalysis.py": ("ffea3ab17ee9c459a0b00f3a44e5e7c23a47e4b01a6e531c417a30c4e0503497",),
    "analysis/__init__.py": ("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",),
    "busco_tools/__init__.py": ("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",),
    "busco_tools/bbtools.py": ("d2674101f0ce78190f96f05976e32b472e6e357a8e7b387b6c68bfbf08ef9ec4", "9a5659eeac4a28b95c1f87c883e96439f740300a6db98b1dd18a191eeb1b7b02"),
    "run_BUSCO.py": ("0193057c468d8e49a1071b24134cd3342cdc5bf82a340436e1953e408e7dba9c", "52d2f3da2467f75440cc8edf36d87864a165dfae98d21aaa126b93b00f4332bb"),
})
FOOTER = "\n# HybAs retention: module-local I/O; original source is inventoried.\nfrom busco_retention import install as _hybas_install_retention\n_hybas_install_retention(globals())\n"
ENTRY = "if __name__ == '__main__':\n    from busco import run_BUSCO\n    run_BUSCO.main()\n"


class UnsupportedBuscoError(RuntimeError):
    """The installed driver does not match the reviewed I/O contract."""


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_directory():
    spec = importlib.util.find_spec("busco")
    if spec is None or not spec.submodule_search_locations:
        raise UnsupportedBuscoError("BUSCO 6.0.0 is not installed in this interpreter")
    return Path(next(iter(spec.submodule_search_locations))).resolve()


def _verify_bytes(name, content):
    expected = SOURCE_HASHES.get(name, ())
    if hashlib.sha256(content).hexdigest() in expected:
        return
    if name == "ConfigManager.py":
        # Conda substitutes this one literal prefix during installation. Check
        # the whole remainder against the reviewed patch; never rewrite it.
        normalized, count = re.subn(
            rb'(?m)^            or "[^"\n]+/share/busco/config\.ini"$',
            b'            or "/opt/anaconda1anaconda2anaconda3/share/busco/config.ini"',
            content,
        )
        if count == 1 and hashlib.sha256(normalized).hexdigest() in expected:
            return
    raise UnsupportedBuscoError(f"Unreviewed BUSCO source: {name}")


def verify_source(source):
    """Verify without importing BUSCO, whose imports create process logs."""
    source = Path(source).resolve()
    version = None
    for node in ast.parse((source / "_version.py").read_text()).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__version__"
            for target in node.targets
        ):
            version = ast.literal_eval(node.value)
    if version != BUSCO_VERSION:
        raise UnsupportedBuscoError(f"Expected BUSCO {BUSCO_VERSION}; found {version!r}")
    found = {path.relative_to(source).as_posix() for path in source.rglob("*.py")}
    if found != set(SOURCE_HASHES):
        raise UnsupportedBuscoError("BUSCO package has missing or additional unreviewed Python modules")
    for name in SOURCE_HASHES:
        path = source / name
        if path.is_symlink() or not path.is_file():
            raise UnsupportedBuscoError(f"Unreviewed BUSCO source: {name}")
        _verify_bytes(name, path.read_bytes())
    return source


class CapturedSubprocess:
    """Keep native calls and returned/raised objects unchanged; retain streams."""

    def __init__(self, native, session):
        self._native, self._session = native, session

    def __getattr__(self, name):
        return getattr(self._native, name)

    def run(self, *args, **kwargs):
        command = args[0] if args else kwargs.get("args")
        try:
            result = self._native.run(*args, **kwargs)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            self._session.capture_process(command, error)
            raise
        self._session.capture_process(command, result)
        return result


class JobLogOpen:
    """Preserve legacy log baselines once, then retain each complete job stream.

    BUSCO deliberately shares these wb paths between parallel jobs. Snapshotting
    before every truncate races with another job's writes. A short cross-process
    lock protects the first baseline snapshot/open; commands still run in
    parallel, and CapturedSubprocess saves every stream before native code writes
    it into the legacy shared file. No native command or file handle is changed.
    """

    def __init__(self, session):
        self.session = session

    def __call__(self, file, mode="r", *args, **kwargs):
        caller = sys._getframe(1)
        job = caller.f_locals.get("self")
        if (mode != "wb" or caller.f_code.co_name != "run"
                or caller.f_globals.get("__name__") != "busco.busco_tools.Toolset"
                or file not in (getattr(job, "job_outlogger", None), getattr(job, "job_errlogger", None))):
            return self.session.open(file, mode, *args, **kwargs)
        path = self.session._path(file, writable=True)
        key = hashlib.sha256(os.fsencode(path)).hexdigest()
        marker = Path(self.session.path) / "busco-log-opens" / key
        try:
            fd = os.open(marker, os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0), 0o600)
            with os.fdopen(fd, "r+b") as lock:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
                state = lock.read()
                if state not in (b"", b"1"):
                    raise RuntimeError("Invalid BUSCO legacy-log retention state")
                if state == b"1":
                    return builtins.open(file, mode, *args, **kwargs)
                handle = self.session.open(file, mode, *args, **kwargs)
                try:
                    lock.write(b"1")
                    lock.flush()
                    os.fsync(lock.fileno())
                except BaseException:
                    handle.close()
                    raise
                return handle
        except Exception as error:
            self.session._fail("BUSCO legacy-log preservation failed: " + str(path), error)


class ReviewedOS:
    """Keep the known BUSCO absent-archive cleanup behavior without deletion."""

    def __init__(self, session, module_name):
        self._proxy = session.os_proxy()
        self._module_name = module_name

    def __getattr__(self, name):
        return getattr(self._proxy, name)

    def remove(self, path, *args, **kwargs):
        caller = sys._getframe(1)
        if (self._module_name == "busco.analysis.BuscoAnalysis"
                and caller.f_code.co_name == "cleanup"
                and caller.f_globals.get("__name__") == self._module_name
                and not args and not kwargs and not os.path.lexists(path)):
            # Upstream passes a bare database-archive filename here, relative to
            # cwd. Do not turn its harmless missing-file case into a failure,
            # or call remove after the check (which would introduce a race).
            raise FileNotFoundError(2, os.strerror(2), os.fspath(path))
        return self._proxy.remove(path, *args, **kwargs)


class ReviewedShutil:
    """Allow two reviewed moves of generated external files into run output."""

    def __init__(self, session, module_name):
        self._session = session
        self._proxy = session.shutil_proxy()
        self._module_name = module_name

    def __getattr__(self, name):
        return getattr(self._proxy, name)

    def move(self, *args, **kwargs):
        caller = sys._getframe(1)
        reviewed = {
            ("busco.BuscoRunner", "move_log_file"),
            ("busco.busco_tools.augustus", "move_retraining_parameters"),
        }
        if (self._module_name, caller.f_code.co_name) in reviewed and caller.f_globals.get("__name__") == self._module_name:
            return self._session.move_preserving_destination(*args, **kwargs)
        return self._proxy.move(*args, **kwargs)


def install(namespace):
    """Called by reviewed private modules, including under multiprocessing spawn."""
    if namespace.get("_hybas_retention_installed"):
        return
    session_path = os.environ.get("HYBAS_NATIVE_RETENTION_SESSION")
    if not session_path:
        raise RuntimeError("Private BUSCO driver requires its retention session")
    session = RetentionSession.attach(session_path)
    if "os" in namespace:
        namespace["os"] = ReviewedOS(session, namespace.get("__name__"))
    if "shutil" in namespace:
        namespace["shutil"] = ReviewedShutil(session, namespace.get("__name__"))
    namespace["open"] = session.open
    if namespace.get("__name__") == "busco.busco_tools.Toolset":
        namespace["subprocess"] = CapturedSubprocess(namespace["subprocess"], session)
        namespace["open"] = JobLogOpen(session)
    namespace["_hybas_retention_installed"] = True


def prepare(source, session, arguments):
    """Create a private driver copy; never patch a shared installed package."""
    source = verify_source(source)
    control = Path(session.path)
    (control / "busco-log-opens").mkdir()
    package_parent = control / "driver"
    package = package_parent / "busco"
    package.mkdir(parents=True)
    inventory = {}
    for item in sorted(source.rglob("*")):
        relative = item.relative_to(source)
        if "__pycache__" in relative.parts or item.suffix in {".pyc", ".pyo"}:
            continue
        if item.is_symlink():
            raise UnsupportedBuscoError(f"BUSCO package contains a symlink: {relative}")
        if not item.is_file():
            continue
        destination = package / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        original = item.read_bytes()
        digest = hashlib.sha256(original).hexdigest()
        if item.suffix == ".py":
            _verify_bytes(relative.as_posix(), original)
        elif item.suffix in {".so", ".pyd", ".dll"}:
            raise UnsupportedBuscoError(f"Unreviewed BUSCO extension module: {relative}")
        reviewed = REVIEWED_SOURCES.get(relative.as_posix())
        if reviewed is not None and digest != reviewed:
            raise UnsupportedBuscoError(f"BUSCO source changed during copy: {relative}")
        copied = original + FOOTER.encode() if reviewed is not None else original
        destination.write_bytes(copied)
        inventory[relative.as_posix()] = {
            "source_sha256": digest,
            "private_sha256": hashlib.sha256(copied).hexdigest(),
            "retention_adapted": reviewed is not None,
        }
    support = control / "support"
    support.mkdir()
    wrapper_inventory = {}
    for name in ("busco_retention.py", "native_output_retention.py"):
        destination = support / name
        shutil.copy2(Path(__file__).with_name(name), destination)
        wrapper_inventory[name] = sha256(destination)
    entry = control / "busco-entry.py"
    entry.write_text(ENTRY)
    receipt = {
        "schema": "hybas-busco-retention-v1",
        "busco_version": BUSCO_VERSION,
        "reviewed_upstream_commit": UPSTREAM_COMMIT,
        "installed_source": str(source),
        "arguments": list(arguments),
        "source_files": inventory,
        "adapter_files": wrapper_inventory,
        "entry_sha256": sha256(entry),
        "status": "prepared",
        "scope": "Reviewed BUSCO Python deletion/replacement I/O and native job stdout/stderr; native subprocess internal scratch is not intercepted",
    }
    receipt_path = control / "busco-driver.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    env = os.environ.copy()
    env.update(session.env)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_parent), str(support)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
    )
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return entry, env, receipt_path, receipt


def _option(arguments, *names):
    value = None
    for index, item in enumerate(arguments):
        for name in names:
            if item == name and index + 1 < len(arguments):
                value = arguments[index + 1]
            elif item.startswith(name + "="):
                value = item.split("=", 1)[1]
    return value


def launch(arguments, output_root, *, source=None, runner=None):
    """Reject caught preservation failures at the external process boundary."""
    arguments = list(arguments)
    # Preserve lexical identity so the session rejects symlinked ancestors.
    root = Path(os.path.abspath(output_root))
    out_path, out_name = _option(arguments, "--out_path"), _option(arguments, "-o", "--out")
    if out_path is None or not out_name:
        raise ValueError("BUSCO retention requires explicit --out_path and -o/--out")
    main_out = (Path(out_path).resolve() / out_name).resolve()
    if not main_out.is_relative_to(root) or main_out == root or main_out.is_relative_to(root / "native-retained"):
        raise ValueError("BUSCO output must be a separate directory inside --output-root")
    source = verify_source(source if source is not None else source_directory())
    session = RetentionSession.create(root, tool="busco")
    try:
        entry, env, receipt_path, receipt = prepare(source, session, arguments)
        try:
            completed = (runner or subprocess.run)([sys.executable, str(entry), *arguments], env=env)
            receipt["returncode"] = completed.returncode
            receipt["status"] = "completed" if completed.returncode == 0 else "failed"
        except BaseException as error:
            receipt["status"] = "failed"
            receipt["launcher_error"] = str(error)
            raise
        finally:
            try:
                session.validate()
            except BaseException as error:
                receipt["status"] = "failed"
                receipt["preservation_error"] = str(error)
                raise
            finally:
                receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        return completed.returncode
    finally:
        session.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    options = parser.parse_args(argv)
    arguments = options.arguments
    if arguments[:1] == ["--"]:
        arguments = arguments[1:]
    return launch(arguments, options.output_root)


if __name__ == "__main__":
    raise SystemExit(main())
