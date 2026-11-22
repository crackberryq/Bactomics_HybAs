#!/usr/bin/env bash
set -euo pipefail
# HYBAS_BROWSER_PREPARATION_CONTRACT=1
# Run inside Snakemake's environment; never activate another environment here.
for tool in fastqc NanoPlot kraken2 extract_kraken_reads.py busco quast.py seqkit samtools pigz java prodigal dnadiff nucmer; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
java -version >/dev/null 2>&1
samtools --version >/dev/null
python - <<'PY'
from importlib import metadata
import json
import os
from pathlib import Path
import random
import struct
import subprocess
from tempfile import TemporaryDirectory
from choreographer.browsers.chromium import Chromium

# NanoPlot 1.46.2 with Plotly 6.1.2 does not install a browser. Only explicit
# environment installation/repair may prepare a missing browser; checks may not.
expected = {'nanoplot': '1.46.2', 'plotly': '6.1.2', 'kaleido': '1.0.0'}
actual = {name: metadata.version(name) for name in expected}
if actual != expected:
    raise RuntimeError('NanoPlot rendering contract differs: ' + json.dumps(actual))
browser = os.environ.get('BROWSER_PATH') or Chromium.find_browser(skip_local=False)
if not browser and os.environ.get('HYBAS_PREPARE_BROWSER') == '1':
    print('Explicit installation/repair: downloading the required Chrome browser through Kaleido.', flush=True)
    try:
        import kaleido
        browser = str(kaleido.get_chrome_sync())
    except Exception as exc:
        raise RuntimeError('Chrome preparation failed. Check network access and browser dependencies, '
                           'then retry Install/Repair environments. No PNG check has passed.') from exc
if not browser:
    raise RuntimeError('NanoPlot needs Chrome/Chromium for Kaleido 1 PNG export. '
                       'Provision a browser explicitly, set BROWSER_PATH if needed, then rerun this check. '
                       'Use an explicit Install/Repair action, or provision Chrome/Chromium yourself. '
                       'No browser is downloaded during checks.')
browser_version = subprocess.check_output([browser, '--version'], text=True, timeout=15).strip()
print('NanoPlot rendering runtime: ' + json.dumps({'packages': actual, 'browser_path': browser,
      'browser_version': browser_version}, sort_keys=True), flush=True)
with TemporaryDirectory(prefix='hybas-nanoplot-check-') as folder:
    root = Path(folder)
    rng = random.Random(42)
    # Synthetic software fixture: varying lengths/quality exercise dot AND KDE.
    # Uniform quality or length would skip the faulty path and falsely pass.
    with (root / 'synthetic.fastq').open('w') as handle:
        for index in range(120):
            length, quality = rng.randrange(300, 1600), rng.randrange(8, 30)
            sequence = ('ACGT' * ((length + 3) // 4))[:length]
            handle.write(f'@synthetic_{index}\n{sequence}\n+\n' + chr(33 + quality) * length + '\n')
    env = {**os.environ, 'BROWSER_PATH': str(browser)}
    subprocess.run(['NanoPlot', '--fastq', str(root / 'synthetic.fastq'), '--outdir', str(root / 'out'),
                    '--threads', '1', '--plots', 'dot', 'kde', '--dpi', '72', '--tsv_stats', '--prefix', 'health_'],
                   check=True, timeout=180, env=env)
    report = root / 'out/health_NanoPlot-report.html'
    if not report.is_file() or 'Plotly.newPlot' not in report.read_text():
        raise RuntimeError('NanoPlot did not generate its HTML report')
    for name in ('health_LengthvsQualityScatterPlot_dot', 'health_LengthvsQualityScatterPlot_kde'):
        html, png = root / 'out' / (name + '.html'), root / 'out' / (name + '.png')
        if not html.is_file() or 'Plotly.newPlot' not in html.read_text():
            raise RuntimeError('NanoPlot did not generate ' + name + '.html')
        data = png.read_bytes() if png.is_file() else b''
        if (len(data) < 45 or not data.startswith(b'\x89PNG\r\n\x1a\n')
                or data[12:16] != b'IHDR' or not all(struct.unpack('>II', data[16:24]))
                or data[-12:] != b'\x00\x00\x00\x00IEND\xaeB`\x82'):
            raise RuntimeError('NanoPlot did not export a complete PNG for ' + name + '. '
                               'Check the browser/renderer log; HTML alone is not a passing check.')
print('PASS: NanoPlot synthetic FASTQ, HTML report, dot/KDE HTML and complete PNG exports')
PY
