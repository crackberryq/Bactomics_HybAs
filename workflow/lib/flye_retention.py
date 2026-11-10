"""Preserve reviewed Flye 2.9.6 products without changing its analysis.

Only module-local filesystem bindings are replaced. Native cleanup still makes
its original paths disappear, and native log/state writes still replace their
active contents. Earlier bytes are retained by RetentionSession first.
Native subprocess scratch (including samtools sort spills) is outside this
adapter's scope. Each invocation requires a new, empty output directory.
"""
from __future__ import annotations

import argparse
import builtins
from contextlib import contextmanager
import hashlib
import importlib
import os
from pathlib import Path
import sys

from native_output_retention import RetentionError, RetentionSession


SUPPORTED_VERSION = "2.9.6"
UPSTREAM_COMMIT = "886b8c17412cdf3a2868a28237bca6c5ad1da156"
MODULE_SHA256 = {
    "flye.main": "9fa37d47e9b7c4c0d9543c94aa6c59eb2a564c11bc36e54aab5feb22f4a221ab",
    "flye.polishing.polish": "c1ab3569740f19872a9a8fe2b7bc182ffae7534772d704486270e110f6cfb27e",
    "flye.polishing.alignment": "3862efb29d4e64d89b10f01c935d19e763eb15de53580506f6ec7f784254e4e7",
}
ADAPTER_SCHEMA = "hybas-flye-native-output-retention-v1"


def verified_modules():
    """Refuse an unreviewed installation before any native workflow executes."""
    version = importlib.import_module("flye.__version__").__version__
    if version != SUPPORTED_VERSION:
        raise RetentionError(
            f"Flye retention requires reviewed version {SUPPORTED_VERSION}; found {version!r}"
        )
    modules = {}
    for name, expected in MODULE_SHA256.items():
        module = importlib.import_module(name)
        source = Path(module.__file__)
        if source.suffix != ".py" or hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise RetentionError(f"Flye retention source guard failed for {name}: {source}")
        modules[name] = module
    if modules["flye.main"].__version__ != version:
        raise RetentionError("Conflicting Flye runtime version identities")
    return modules


def _fresh_root(root, native_args):
    root = Path(os.path.abspath(os.fspath(root)))
    if root.resolve() != root or not root.is_dir() or any(root.iterdir()):
        raise RetentionError("Flye retention requires an existing, empty output directory without symlinked ancestors")
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("-o", "--out-dir", required=True)
    args, _ = parser.parse_known_args(native_args)
    native_root = Path(os.path.abspath(args.out_dir))
    if native_root != root:
        raise RetentionError("Flye --out-dir must match the retention --root")
    return root


def _reviewed_open(session, basename):
    def retained_open(file, mode="r", *args, **kwargs):
        # Native /dev/null and ordinary reads retain normal Python semantics.
        # Destructive opens only occur at the reviewed, hash-bound call sites.
        if isinstance(file, (str, bytes, os.PathLike)) and "w" in mode:
            if os.path.basename(os.fsdecode(file)) == basename:
                return session.open(file, mode, *args, **kwargs)
        return builtins.open(file, mode, *args, **kwargs)
    return retained_open


@contextmanager
def preserve_native_outputs(modules, session):
    """Install scoped adapters; restore bindings even after a native failure."""
    replacements = [
        (modules["flye.main"], "os", session.os_proxy()),
        (modules["flye.polishing.polish"], "os", session.os_proxy()),
        (modules["flye.main"], "open", _reviewed_open(session, "params.json")),
        (modules["flye.polishing.alignment"], "open", _reviewed_open(session, "minimap.stderr")),
    ]
    absent = object()
    originals = []
    try:
        for module, name, replacement in replacements:
            originals.append((module, name, module.__dict__.get(name, absent)))
            setattr(module, name, replacement)
        yield
    finally:
        for module, name, original in reversed(originals):
            if original is absent:
                delattr(module, name)
            else:
                setattr(module, name, original)
        # An upstream catch must never turn failed preservation into success.
        session.validate()


def run(root, native_args):
    """Run the native entry point with the original argument list unchanged."""
    native_args = list(native_args)
    root = _fresh_root(root, native_args)
    modules = verified_modules()
    metadata = {
        "adapter_schema": ADAPTER_SCHEMA,
        "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "native_version": SUPPORTED_VERSION,
        "native_source_commit": UPSTREAM_COMMIT,
        "native_module_sha256": dict(MODULE_SHA256),
        "native_argv": ["flye", *native_args],
        "scope": "reviewed Python cleanup; mapping stderr and params.json revisions; native subprocess scratch excluded",
    }
    session = RetentionSession.create(root, tool="flye", metadata=metadata)
    original_argv = sys.argv
    try:
        with preserve_native_outputs(modules, session):
            sys.argv = ["flye", *native_args]
            return modules["flye.main"].main()
    finally:
        sys.argv = original_argv
        session.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="Fresh Flye output directory")
    parser.add_argument("native_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if not args.native_args or args.native_args[0] != "--":
        parser.error("Separate the unchanged Flye arguments with --")
    try:
        return run(args.root, args.native_args[1:])
    except (RetentionError, OSError, ImportError) as error:
        print(f"Flye output retention failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
