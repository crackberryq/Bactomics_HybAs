#!/usr/bin/env python3
"""Export an experimental read-preservation bundle from a HybAs review receipt."""
import argparse
import json
from pathlib import Path
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lib'))
from aware_plus import preserve

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', required=True, help='HybAs identification review.json')
    parser.add_argument('--output', required=True, help='New output directory; existing results are never replaced')
    args = parser.parse_args()
    started, cpu = time.monotonic(), time.process_time()
    status = 'failed'
    try:
        result = preserve(args.review, args.output)
        status = 'succeeded'
        print(result['version'] + ': preserved ' + str(result['output_units']))
    finally:
        # The command log retains metrics even when no output bundle can be published.
        print(json.dumps({'event': 'aware_plus_preservation_metrics', 'status': status,
            'wall_seconds': time.monotonic()-started, 'cpu_seconds': time.process_time()-cpu,
            'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == 'darwin' else 1)}), file=sys.stderr)
