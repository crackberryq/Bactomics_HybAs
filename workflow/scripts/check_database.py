"""Validate database structure before committing a native download rule."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from workflow_support import busco_db_ready

parser = argparse.ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()
if not busco_db_ready(args.path):
    parser.exit(1, f"Incomplete BUSCO dataset at {args.path}: dataset.cfg, hmms, scores_cutoff and lengths_cutoff are required\n")
