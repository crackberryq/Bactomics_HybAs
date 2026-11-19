"""Verify a BUSCO lineage or prepare it in the workflow's stable location."""
import argparse
import fcntl
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from workflow_support import busco_db_ready


def publish(source, target):
    # Separate admin workspaces can prepare the same shared lineage. Serialize
    # only the final directory change and reuse another completed download.
    fd = os.open(target.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        if busco_db_ready(target):
            return
        backup = None
        if target.is_symlink() or (target.exists() and (not target.is_dir() or any(target.iterdir()))):
            saved = Path(tempfile.mkdtemp(prefix=f".previous-incomplete.{target.name}.", dir=target.parent))
            backup = saved / target.name
            target.rename(backup)
            print(f"Preserved previous incomplete dataset at {backup}", flush=True)
        try:
            # The staging directory uses the same filesystem as the target.
            source.replace(target)
        except OSError:
            if backup is not None:
                backup.rename(target)
            raise
    finally:
        os.close(fd)


def prepare(target):
    if busco_db_ready(target):
        print(f"Verified existing BUSCO dataset at {target}")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".download.{target.name}.", dir=target.parent) as tmp:
        root = Path(tmp)
        subprocess.run(["busco", "--download", target.name, "--download_path", str(root)], check=True)
        # BUSCO's downloader groups datasets by type. Also recognize the flat
        # layout accepted by earlier HybAs versions.
        candidates = [root / "lineages" / target.name, root / target.name]
        for candidate in candidates:
            if busco_db_ready(candidate):
                publish(candidate, target)
                print(f"Verified BUSCO dataset at {target}")
                return
        paths = " and ".join(map(str, candidates))
        raise ValueError(f"Incomplete BUSCO download: checked {paths}. "
                         "Expected dataset.cfg, hmms/*.hmm, scores_cutoff and lengths_cutoff.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--prepare", action="store_true", help="Download and verify the lineage if needed")
    args = parser.parse_args()
    if args.prepare:
        try:
            prepare(args.path)
        except (OSError, ValueError, subprocess.CalledProcessError) as exc:
            parser.exit(1, f"{exc}\n")
    elif not busco_db_ready(args.path):
        parser.exit(1, f"Incomplete BUSCO dataset at {args.path}: dataset.cfg, hmms, scores_cutoff and lengths_cutoff are required\n")
