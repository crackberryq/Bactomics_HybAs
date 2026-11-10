"""Retain local products before Snakemake prepares a replacement job.

The engine removes output files and directory() trees *before* the rule body.
This small, version-gated adapter runs before that removal, including retries.
Only the installed workflow is affected. No analysis parameters are changed.
"""
from __future__ import annotations

from datetime import datetime, timezone
from functools import wraps
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import re
import stat
import uuid


SUPPORTED_SNAKEMAKE = {"9.6.2", "9.27.1"}


def _write_manifest(path, value):
    temporary = path.with_suffix(".pending")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def _absolute(path):
    # Do not follow the final symlink: retain the link, never move its target.
    return Path(os.path.abspath(os.fspath(path)))


def _reject_hardlinks(path):
    """A rename cannot isolate bytes that another hardlink can still change."""
    pending = [path]
    while pending:
        current = pending.pop()
        info = current.lstat()  # Never follow a symbolic link or traverse its target.
        if stat.S_ISREG(info.st_mode) and info.st_nlink > 1:
            raise ValueError(f"Retention refuses a hardlinked product: {current}")
        if stat.S_ISDIR(info.st_mode):
            with os.scandir(current) as entries:
                pending.extend(Path(entry.path) for entry in entries)


def _validated_paths(base, paths, inputs=()):
    """Validate the entire owned selection before moving any of its products."""
    archive_root = base / "retained"
    candidates = []
    for value in paths:
        path = _absolute(value)
        if path == base or not path.is_relative_to(base) or path.is_relative_to(archive_root):
            raise ValueError(f"Refusing to retain an unowned output path: {path}")
        # A symlinked ancestor could turn a lexical child into external data.
        if path.parent.resolve() != path.parent or base.resolve() != base:
            raise ValueError(f"Retention requires output directories without symlinked ancestors: {path}")
        for value_in in inputs:
            source = _absolute(value_in)
            if source == path or source.is_relative_to(path):
                raise ValueError(f"Refusing to move a product containing a current input: {path}")
        if os.path.lexists(path):
            candidates.append(path)
    # Retaining a directory already retains all of its listed children.
    selected = []
    for path in sorted(set(candidates), key=lambda item: (len(item.parts), str(item))):
        if not any(path.is_relative_to(parent) for parent in selected):
            selected.append(path)
    if archive_root.is_symlink():
        raise ValueError(f"Retention archive must not be a symlink: {archive_root}")
    for path in selected:
        _reject_hardlinks(path)
    return selected


def _archive_selected(base, selected, *, rule, wildcards=None, reason="before-rerun"):
    if not selected:
        return None
    archive_root = base / "retained"
    label = re.sub(r"[^A-Za-z0-9_.-]", "_", str(rule))[:80]
    identifier = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + "-" + label + "-" + uuid.uuid4().hex
    destination = archive_root / identifier
    destination.mkdir(parents=True, exist_ok=False)
    manifest = {"schema": "workflow-output-retention-v1", "state": "archiving",
        "reason": reason, "rule": str(rule), "wildcards": dict(wildcards or {}),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "prior_execution_status": "not_inferred", "entries": []}
    for path in selected:
        relative = path.relative_to(base).as_posix()
        entry = {"original": str(path), "archive": "files/" + relative, "moved": False}
        if path.is_symlink():
            entry.update(kind="symlink", target=os.readlink(path))
        else:
            entry["kind"] = "directory" if path.is_dir() else "file"
        manifest["entries"].append(entry)
    manifest_path = destination / "manifest.json"
    _write_manifest(manifest_path, manifest)
    try:
        for path, entry in zip(selected, manifest["entries"]):
            target = destination / entry["archive"]
            target.parent.mkdir(parents=True, exist_ok=True)
            path.rename(target)
            entry["moved"] = True
            _write_manifest(manifest_path, manifest)
        manifest["state"] = "complete"
        _write_manifest(manifest_path, manifest)
    except BaseException as error:
        manifest.update(state="interrupted", error=str(error))
        _write_manifest(manifest_path, manifest)
        raise RuntimeError(f"Output preservation failed; no replacement job may run. Inspect {manifest_path}") from error
    return destination


def archive_paths(base_dir, paths, *, rule, wildcards=None, inputs=(), reason="before-rerun"):
    """Move owned local products to a unique archive; fail before any deletion.

    Moves preserve bytes without hardlinks or a second full-size copy. Existing
    hardlinks are refused because a rename would not isolate their bytes. The
    manifest is written before the first move and updated after every move, so
    interrupted preservation remains inspectable. Cross-filesystem moves fail
    closed; there is deliberately no copy/delete fallback.
    """
    base = _absolute(base_dir)
    selected = _validated_paths(base, paths, inputs)
    return _archive_selected(base, selected, rule=rule, wildcards=wildcards, reason=reason)


def _archive_job_paths(base_dir, owned_roots, paths, *, rule, wildcards=None, inputs=()):
    """Validate every root's selection first; archive each under its owning root."""
    roots = sorted({_absolute(base_dir), *map(_absolute, owned_roots)},
                   key=lambda root: (-len(root.parts), str(root)))
    if any(root.is_relative_to(other / "retained") for root in roots for other in roots if root != other):
        raise ValueError("An owned root cannot be inside another root's retention archive")
    groups = {root: [] for root in roots}
    for value in paths:
        path = _absolute(value)
        owner = next((root for root in roots if path.is_relative_to(root)), None)
        if owner is None:
            raise ValueError(f"Refusing to retain an unowned output path: {path}")
        # A broad directory must not swallow a separately registered root or
        # its retained history. Each registered generated root is independent.
        if any(root != owner and root.is_relative_to(path) for root in roots):
            raise ValueError(f"Refusing to move a product containing another owned root: {path}")
        groups[owner].append(path)
    validated = [(root, _validated_paths(root, group, inputs)) for root, group in groups.items()]
    return [_archive_selected(root, selected, rule=rule, wildcards=wildcards)
            for root, selected in validated if selected]


def _job_paths(job):
    products = list(job.output) + list(job.log)
    if job.benchmark:
        products.append(job.benchmark)
    # Native tool sidecars are not necessarily individual Snakemake outputs.
    products.extend(getattr(job.params, "retention_paths", []))
    return products


def install(workflow, base_dir, owned_roots=()):
    """Install before scheduling, including native standalone invocation.

    ``params.retention_paths`` lists directories/files exclusively owned by a
    rule. Do not list a shared work directory or a rule's inputs there.
    ``owned_roots`` explicitly registers additional generated roots, such as a
    database root. Their archives stay in that root's ``retained/`` directory.
    """
    from snakemake.jobs import Job
    version = importlib.metadata.version("snakemake")
    if version not in SUPPORTED_SNAKEMAKE:
        raise RuntimeError(f"Output retention has not been verified for Snakemake {version}; use one of {sorted(SUPPORTED_SNAKEMAKE)}")
    original = Job.remove_existing_output
    if not inspect.iscoroutinefunction(original) or list(inspect.signature(original).parameters) != ["self"]:
        raise RuntimeError("Snakemake output preparation changed; refusing to run without output retention")
    if workflow.execution_settings is not None:
        workflow.execution_settings.keep_incomplete = True
    if workflow.storage_settings is not None:
        workflow.storage_settings.notemp = True
    workflow._artifact_retention_base = str(_absolute(base_dir))
    workflow._artifact_retention_roots = tuple(str(_absolute(root)) for root in owned_roots)
    if getattr(original, "_artifact_retention", False):
        return

    @wraps(original)
    async def preserve_then_remove(self):
        owner = self.rule.workflow
        base = getattr(owner, "_artifact_retention_base", None)
        if base is not None:
            owner.execution_settings.keep_incomplete = True
            owner.storage_settings.notemp = True
            # Fail closed rather than accidentally delete a remote or in-place
            # product; shipped workflows use ordinary local output paths.
            from snakemake.io import is_flagged
            if any(getattr(item, "is_storage", False) or is_flagged(item, "update") for item in self.output):
                raise RuntimeError("Output retention requires local, non-in-place products")
            _archive_job_paths(base, getattr(owner, "_artifact_retention_roots", ()),
                               _job_paths(self), rule=self.rule.name,
                               wildcards=self.wildcards_dict, inputs=self.input)
        await original(self)

    preserve_then_remove._artifact_retention = True
    Job.remove_existing_output = preserve_then_remove
