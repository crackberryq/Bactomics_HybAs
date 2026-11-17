# QC rendering environments

NanoPlot 1.46.2 declares `kaleido>=1.0.0`; MultiQC 1.32 declares
`kaleido==0.2.1`. They cannot share a dependency-consistent environment.
The QC recipe keeps NanoPlot with Plotly 6.1.2 and Kaleido 1.0.0; the separate
MultiQC recipe uses Plotly 6.0.1 and Kaleido 0.2.1. Plotly 7 removed the
`create_2d_density` API used by NanoPlot 1.46.2. These are software dependency
constraints, not evidence of analytical validity.

Kaleido 1 needs a working Chrome or Chromium installation. A supported
Bactologix **Install** or **Repair** action can prepare a missing browser through
the installed Kaleido package. The post-deployment hook permits that download
only when explicitly enabled with `HYBAS_PREPARE_BROWSER=1`. Existing browsers
are reused. A **Check** never enables preparation. Native CLI users can opt in
with the same flag when creating environments, or provision a browser through
their system package manager or `kaleido_get_chrome`. Downloading a browser
requires network access. Set `BROWSER_PATH`
to an existing browser executable when automatic discovery is unsuitable.
The QC environment manifest includes a `qc.graphics_runtime.json` observation
with the discovered browser path, executable SHA256 and version output. The
MultiQC manifest records its bundled Kaleido renderer separately. These
observations do not establish that a plot was produced, and hashing an executable
does not fingerprint all shared libraries or the target of a shell wrapper.
Environment checks do not install a browser or silently omit PNG exports.

The post-deployment checks execute real NanoPlot on deterministic synthetic
FASTQ reads with varying length and quality, and real MultiQC on synthetic
custom content. Both must produce the expected HTML and complete PNG files;
an executable's version output is insufficient. The native-environment CI
job explicitly enables browser preparation in the native NanoPlot hook and runs both checks
in fresh Conda environments. A source change is not qualification until those
clean environment checks pass on the intended platform.

Use a new pipeline installation/cache and a revised draft for the split
recipes. Keep original run configuration, completion snapshots, outputs and
environment manifests intact. Historical repairs require separately dated
before/after package records; a current package version does not establish
which version produced an older result.
