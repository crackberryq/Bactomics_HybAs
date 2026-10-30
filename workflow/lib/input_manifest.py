"""HybAs v10 metadata contract. No sequence processing or workflow execution.

This deliberately has no dependency on the v9 controller. A valid manifest is
not an executable plan or evidence that a proposed route has been qualified.
"""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re


SCHEMA_VERSION = 1
MODES = (
    "illumina_only", "dual_illumina_runs", "dual_illumina_complementary",
    "ont_only", "illumina_ont", "pacbio_only", "illumina_pacbio",
)
FASTQ_SUFFIXES = (".fq", ".fastq", ".fq.gz", ".fastq.gz")


class ManifestError(ValueError):
    """Invalid or ambiguous input metadata, with a location for the user."""


def _error(where, message):
    raise ManifestError(f"{where}: {message}")


def _object(value, required, optional, where):
    if not isinstance(value, dict) or any(not isinstance(k, str) for k in value):
        _error(where, "expected an object with string keys")
    missing = set(required) - value.keys()
    unknown = value.keys() - set(required) - set(optional)
    if missing:
        _error(where, f"missing fields: {', '.join(sorted(missing))}")
    if unknown:
        _error(where, f"unknown fields: {', '.join(sorted(unknown))}")


def _items(value, where):
    if not isinstance(value, list) or not value:
        _error(where, "expected a non-empty list")
    return value


def _choice(value, choices, where):
    if not isinstance(value, str) or value not in choices:
        _error(where, f"expected one of: {', '.join(choices)}")
    return value


def _identifier(value, where):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        _error(where, "use an identifier starting with a letter or digit; then letters, digits, '.', '_' or '-'")
    return value


def _unique(value, seen, where):
    if value in seen:
        _error(where, f"duplicate identifier {value!r}")
    seen.add(value)


def _path(value, base_dir, seen, where):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        _error(where, "expected a non-empty path without surrounding whitespace")
    if any(ord(c) < 32 for c in value) or any(c in value for c in "*?[]"):
        _error(where, "use a literal file path, not a glob or control characters")
    if "://" in value or "\\" in value or re.match(r"^[A-Za-z]:", value):
        _error(where, "use a local Linux/WSL path, not a URL or Windows path")
    if not value.endswith(FASTQ_SUFFIXES):
        _error(where, "supported declarations are .fastq/.fq, optionally .gz")
    try:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = base_dir / path
        path = str(path.resolve())
    except (OSError, RuntimeError, ValueError) as exc:
        _error(where, f"cannot resolve input path: {exc}")
    if path in seen:
        _error(where, f"file is already declared at {seen[path]}")
    seen[path] = where
    return path


def _library(value, base_dir, seen_files, where):
    _object(value, ("library_id", "platform", "read_groups"),
            ("library_type", "orientation", "insert_size_bp", "chemistry", "preprocessing"), where)
    out = deepcopy(value)
    _identifier(out["library_id"], f"{where}.library_id")
    platform = _choice(out["platform"], ("illumina", "ont", "pacbio"), f"{where}.platform")
    if platform == "illumina":
        _choice(out.get("library_type"), ("paired_end", "mate_pair"), f"{where}.library_type")
        _choice(out.get("orientation"), ("fr", "rf", "ff"), f"{where}.orientation")
        if "chemistry" in out:
            _error(where, "chemistry is only used for PacBio in schema v1")
        if "preprocessing" in out:
            _choice(out["preprocessing"], ("fastp", "pretrimmed", "nextera_mate_pair"), f"{where}.preprocessing")
            if out["preprocessing"] == "nextera_mate_pair" and (out["library_type"] != "mate_pair" or out["orientation"] != "rf"):
                _error(where, "Nextera mate-pair preprocessing requires mate_pair/rf input")
        if "insert_size_bp" in out and (type(out["insert_size_bp"]) is not int or out["insert_size_bp"] <= 0):
            _error(where, "insert_size_bp must be a positive integer")
    else:
        if any(key in out for key in ("library_type", "orientation", "insert_size_bp", "preprocessing")):
            _error(where, "paired-library fields are only valid for Illumina")
        if platform == "pacbio":
            _choice(out.get("chemistry"), ("hifi", "clr"), f"{where}.chemistry")
        elif "chemistry" in out:
            _error(where, "chemistry is only used for PacBio in schema v1")

    groups, seen_groups = [], set()
    for index, group in enumerate(_items(out["read_groups"], f"{where}.read_groups")):
        loc = f"{where}.read_groups[{index}]"
        file_fields = ("r1", "r2") if platform == "illumina" else ("reads",)
        _object(group, ("read_group_id", "run_id", *file_fields), ("lane_id",), loc)
        _unique(_identifier(group["read_group_id"], loc + ".read_group_id"), seen_groups, loc)
        _identifier(group["run_id"], loc + ".run_id")
        if "lane_id" in group:
            _identifier(group["lane_id"], loc + ".lane_id")
        normalized = dict(group)
        for key in file_fields:
            normalized[key] = _path(group[key], base_dir, seen_files, loc + "." + key)
        groups.append(normalized)
    out["read_groups"] = groups
    return out


def _select_mode(sample, where):
    libraries = sample["libraries"]
    platforms = {lib["platform"] for lib in libraries}
    illumina = [lib for lib in libraries if lib["platform"] == "illumina"]
    has_mp = any(lib["library_type"] == "mate_pair" for lib in illumina)
    has_pe = any(lib["library_type"] == "paired_end" for lib in illumina)
    run_ids = {g["run_id"] for lib in illumina for g in lib["read_groups"]}
    pacbio_chemistries = {lib["chemistry"] for lib in libraries if lib["platform"] == "pacbio"}
    if len(pacbio_chemistries) > 1:
        _error(where, "mixed PacBio HiFi and CLR libraries are not supported in schema v1")

    compatible = []
    if platforms == {"illumina"}:
        if has_pe and not has_mp:
            if len(illumina) == 1:
                compatible = ["dual_illumina_runs" if len(run_ids) > 1 else "illumina_only"]
            else:
                compatible = ["dual_illumina_runs", "dual_illumina_complementary"]
        elif has_pe and has_mp:
            compatible = ["dual_illumina_complementary"]
    elif platforms == {"ont"}:
        compatible = ["ont_only"]
    elif platforms == {"pacbio"}:
        compatible = ["pacbio_only"]
    elif platforms == {"illumina", "ont"} and not has_mp:
        compatible = ["illumina_ont"]
    elif platforms == {"illumina", "pacbio"} and not has_mp:
        compatible = ["illumina_pacbio"]
    if not compatible:
        _error(where, "unsupported platform/library combination; no input may be silently ignored")

    requested = sample.get("mode", "auto")
    if requested == "auto":
        if len(compatible) != 1:
            _error(where, "multiple paired-end libraries require an explicit mode: " + ", ".join(compatible))
        selected = compatible[0]
    elif requested not in compatible:
        _error(where, f"mode {requested!r} does not match the declared inputs; compatible: {', '.join(compatible)}")
    else:
        selected = requested

    if selected == "dual_illumina_complementary":
        if any("insert_size_bp" not in lib for lib in illumina):
            _error(where, "complementary libraries require insert_size_bp for every Illumina library")
        if not has_mp and len({lib["insert_size_bp"] for lib in illumina}) < 2:
            _error(where, "complementary paired-end libraries require distinct declared insert sizes")

    return selected, next(iter(pacbio_chemistries), None)


def validate_manifest(document, *, base_dir):
    """Validate declarations and normalize paths, without opening sequence files.

    base_dir is the manifest directory, never the caller's working directory.
    File uniqueness applies across the entire batch after path normalization.
    """
    _object(document, ("schema_version", "samples"), (), "manifest")
    if type(document["schema_version"]) is not int or document["schema_version"] != SCHEMA_VERSION:
        _error("manifest.schema_version", f"expected integer {SCHEMA_VERSION}")
    normalized, selections, seen_samples, seen_files = [], [], set(), {}
    base_dir = Path(base_dir).resolve()
    for index, sample in enumerate(_items(document["samples"], "manifest.samples")):
        where = f"manifest.samples[{index}]"
        _object(sample, ("sample_id", "libraries"), ("mode",), where)
        sample_id = _identifier(sample["sample_id"], where + ".sample_id")
        _unique(sample_id, seen_samples, where)
        requested = _choice(sample.get("mode", "auto"), ("auto", *MODES), where + ".mode")
        libraries, seen_libraries = [], set()
        for lib_index, library in enumerate(_items(sample["libraries"], where + ".libraries")):
            loc = f"{where}.libraries[{lib_index}]"
            result = _library(library, base_dir, seen_files, loc)
            _unique(result["library_id"], seen_libraries, loc)
            libraries.append(result)
        item = {"sample_id": sample_id, "mode": requested, "libraries": libraries}
        selected, chemistry = _select_mode(item, where)
        normalized.append(item)
        selections.append({
            "sample_id": sample_id, "requested_mode": requested, "selected_mode": selected,
            "selection_source": "declared_metadata" if requested == "auto" else "explicit",
            "pacbio_chemistry": chemistry,
            "library_ids": [lib["library_id"] for lib in libraries],
            "execution_ready": False,
            "qualification": "manifest_validation_only",
        })
    return {
        "schema_version": SCHEMA_VERSION, "status": "metadata_validated",
        "execution_ready": False,
        "note": "This preflight validates declarations only. It does not select executable commands, "
                "run tools, verify sequence content, or establish assembly accuracy.",
        "manifest": {"schema_version": SCHEMA_VERSION, "samples": normalized},
        "selections": selections,
    }


def _json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            _error("manifest", f"duplicate key {key!r}")
        result[key] = value
    return result


def load_manifest(path):
    """Read JSON or safe YAML, rejecting duplicate keys rather than hiding them."""
    path = Path(path).resolve()
    raw = path.read_bytes()
    text = raw.decode("utf-8-sig")
    if path.suffix.lower() == ".json":
        try:
            document = json.loads(text, object_pairs_hook=_json_pairs)
        except json.JSONDecodeError as exc:
            raise ManifestError(f"{path}: invalid JSON: {exc}") from exc
    elif path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml
        except ImportError as exc:
            raise ManifestError("YAML manifests require PyYAML; JSON uses the standard library") from exc

        class UniqueSafeLoader(yaml.SafeLoader):
            pass

        def mapping(loader, node):
            loader.flatten_mapping(node)
            pairs = [(loader.construct_object(k), loader.construct_object(v)) for k, v in node.value]
            if any(not isinstance(k, str) for k, _ in pairs):
                _error("manifest", "YAML object keys must be strings")
            return _json_pairs(pairs)

        UniqueSafeLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)
        try:
            document = yaml.load(text, Loader=UniqueSafeLoader)
        except yaml.YAMLError as exc:
            raise ManifestError(f"{path}: invalid YAML: {exc}") from exc
    else:
        raise ManifestError("Manifest filename must end in .json, .yaml or .yml")
    result = validate_manifest(document, base_dir=path.parent)
    result["source"] = {"manifest_path": str(path), "sha256": hashlib.sha256(raw).hexdigest()}
    return result
