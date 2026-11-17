"""Observe rendering dependencies without downloading a browser or creating plots."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import subprocess
import sys


def executable_identity(path):
    selected = Path(path).expanduser().absolute()
    resolved = selected.resolve(strict=True)
    digest = hashlib.sha256()
    with resolved.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(selected), 'resolved_path': str(resolved),
            'executable_sha256': digest.hexdigest(), 'size_bytes': resolved.stat().st_size,
            'scope': 'Selected executable only; shared libraries and wrapper targets are not hashed.'}


def discover_browser():
    configured = os.environ.get('BROWSER_PATH')
    if configured:
        return configured
    # Discovery does not download Chrome. NanoPlot is deliberately not imported.
    from choreographer.browsers.chromium import Chromium
    return Chromium.find_browser(skip_local=False)


def capture(tool):
    names = [tool, 'plotly', 'kaleido', 'choreographer']
    packages = {}
    for name in names:
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None
    receipt = {'schema_version': 1, 'observed_at': datetime.now(timezone.utc).isoformat(),
               'tool': tool, 'python': sys.version, 'packages': packages,
               'scope': 'Current rendering runtime observation; not evidence that any plot was produced.'}
    try:
        if tool == 'nanoplot':
            browser = discover_browser()
            if not browser:
                raise RuntimeError('No installed Chrome/Chromium browser was discovered.')
            identity = executable_identity(browser)
            result = subprocess.run([browser, '--version'], text=True, capture_output=True, timeout=15)
            identity.update(version_output=(result.stdout + result.stderr).strip()[:4096],
                            version_exit_code=result.returncode)
            receipt['browser'] = {'status': 'observed', **identity}
        else:
            path = metadata.distribution('kaleido').locate_file('kaleido/executable/bin/kaleido')
            receipt['bundled_renderer'] = {'status': 'observed', **executable_identity(path)}
    except (OSError, ImportError, RuntimeError, subprocess.SubprocessError) as exc:
        receipt['browser' if tool == 'nanoplot' else 'bundled_renderer'] = {
            'status': 'unavailable', 'reason': str(exc)}
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tool', choices=('nanoplot', 'multiqc'), required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(capture(args.tool), indent=2, sort_keys=True) + '\n')
