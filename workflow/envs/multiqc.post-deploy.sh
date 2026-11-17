#!/usr/bin/env bash
set -euo pipefail
location=$(command -v multiqc) || { echo "Missing executable: multiqc" >&2; exit 1; }
case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: multiqc" >&2; exit 1;; esac
python - <<'PY'
from importlib import metadata
import json
from pathlib import Path
import struct
import subprocess
from tempfile import TemporaryDirectory

expected = {'multiqc': '1.32', 'plotly': '6.0.1', 'kaleido': '0.2.1'}
actual = {name: metadata.version(name) for name in expected}
if actual != expected:
    raise RuntimeError('MultiQC rendering contract differs: ' + json.dumps(actual))
print('MultiQC rendering packages: ' + json.dumps(actual, sort_keys=True), flush=True)
with TemporaryDirectory(prefix='hybas-multiqc-check-') as folder:
    root = Path(folder)
    fixture = {'id': 'hybas_render_check', 'section_name': 'Synthetic rendering check',
               'plot_type': 'bargraph', 'pconfig': {'id': 'hybas_render_check', 'title': 'Synthetic values'},
               'data': {'fixture_a': {'value': 10}, 'fixture_b': {'value': 20}}}
    (root / 'render_mqc.json').write_text(json.dumps(fixture))
    subprocess.run(['multiqc', str(root), '--force', '--export', '--outdir', str(root / 'out')],
                   check=True, timeout=120)
    report = root / 'out/multiqc_report.html'
    if not report.is_file() or 'Synthetic rendering check' not in report.read_text():
        raise RuntimeError('MultiQC did not generate the expected HTML report')
    images = list((root / 'out').rglob('*hybas_render_check*.png'))
    if not images:
        raise RuntimeError('MultiQC did not export a PNG for the synthetic plot')
    for path in images:
        data = path.read_bytes()
        if (len(data) < 45 or not data.startswith(b'\x89PNG\r\n\x1a\n')
                or data[12:16] != b'IHDR' or not all(struct.unpack('>II', data[16:24]))
                or data[-12:] != b'\x00\x00\x00\x00IEND\xaeB`\x82'):
            raise RuntimeError('MultiQC PNG is incomplete: ' + str(path))
print('PASS: MultiQC synthetic HTML report and complete PNG export')
PY
