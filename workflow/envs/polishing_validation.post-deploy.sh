#!/usr/bin/env bash
set -euo pipefail
for tool in bwa samtools dnadiff show-snps nucmer delta-filter; do
    command -v "$tool" >/dev/null
done
python - <<'PY'
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
import pysam

figure = pyplot.figure()
figure.add_subplot(111, projection="3d")
pyplot.close(figure)
print("Polishing validation environment ready")
PY
