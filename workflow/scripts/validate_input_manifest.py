#!/usr/bin/env python3
"""Validate v10 input metadata and print JSON; never launch an assembly."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from input_manifest import ManifestError, load_manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="Versioned JSON or YAML input manifest")
    args = parser.parse_args(argv)
    try:
        result = load_manifest(args.manifest)
    except (ManifestError, OSError, UnicodeError) as exc:
        print(f"Input manifest error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
