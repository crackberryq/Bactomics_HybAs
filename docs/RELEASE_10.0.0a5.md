# HybAs AWARE 10.0.0a5 — Alpha

This is a development candidate until the exact source passes its required
checks and is published as `v10.0.0a5`. A local patch or a passing synthetic
workflow is not a released or scientifically validated pipeline.

## Changes

When Medaka is selected, an early preflight inspects all raw ONT read headers
before merging, filtering or assembly. Automatic selection requires one
resolved basecaller identity throughout the inputs and a matching entry in
the installed Medaka mapping. Missing, mixed, unsupported or mismatched
metadata requires an explicit model decision; unmatched choices require an
acknowledgement, reason and binding to that exact model. No default model is
silently substituted, and an acknowledgement does not establish scientific
suitability. Identification-only endpoints defer polishing review until
continuation.

Each attempt retains its model-resolution command, errors and receipt. Before
consensus, HybAs rechecks the filtered reads, selected model hash and software
identity against preflight. Previous native attempt outputs remain available.
Medaka stays pinned at 2.2.0; this change does not approve an historical model
preset or validate inference accuracy.

The Medaka environment installs its canonical ONT Python dependencies and
pins pysam within Medaka's declared range. PyTorch 2.9.1 and NumPy 2.3.5 stay
Conda-owned; the Python installer is constrained against replacing them.
Setup checks package consistency, namespace ownership, imported native files
and small synthetic library operations. These checks do not load a consensus
model or validate GPU execution. Installed builds and versions remain recorded
for comparison with earlier runs.

NanoPlot 1.46.2 and MultiQC 1.32 now use separate environments because their
Kaleido requirements conflict. The QC environment pins Plotly 6.1.2 and
Kaleido 1.0.0; MultiQC uses Plotly 6.0.1 and Kaleido 0.2.1. The MultiQC
environment is selected only when reporting is enabled, and version records
read its own captured manifest. Analysis thresholds are unchanged, and native
outputs remain retained.

Post-deployment checks run the actual NanoPlot and MultiQC command-line tools
on synthetic inputs and require both HTML and complete PNG output. NanoPlot
checks exercise dot and KDE plots with varying read length and quality.
Chrome or Chromium must be provisioned explicitly for NanoPlot; the checks
do not download a browser. CI provisions Chrome in a separate, visible step.

The QC and MultiQC manifests retain a dated rendering-runtime observation
alongside the Conda and Python records. It includes package versions and the
selected browser or bundled renderer executable identity. Browser discovery
and version inspection do not establish that any plot was produced, and the
receipt states this limit. Existing products, native sidecars, logs,
benchmarks and earlier manifests remain covered by workflow output retention.

## Update and test

Use a new source installation and environment cache, then create a new or
revised draft. Preserve the original source, frozen settings, output files and
completion records when comparing results. Do not repair a historical result
by replacing its recorded package versions with current installations.

Follow [QC environment setup](QC_ENVIRONMENTS.md) to provision the browser
explicitly before checking the QC environment. Missing or unusable rendering
dependencies must be resolved before accepting the environment as ready.
For Bactologix, wait for a compatibility entry reviewed against the exact
published source commit; a new GitHub tag alone does not enable a runner.

## Verification and scope

Publication requires the workflow regressions on Python 3.11 and 3.12,
report-browser checks, native controller and output-retention contracts, a
fresh Medaka environment that passes its real CLI/model-interface checks, and
fresh QC and MultiQC Conda environments that pass real HTML/PNG rendering
checks. Medaka environment checks download no model weights and run no
inference. The release notes append the tested source commit and CI run URL at
publication. Published tags retain their original source identities.

Synthetic workflow tests exercise scheduling, output retention and provenance
with tool stubs; they do not run a complete biological analysis. Rendering
checks exercise the actual CLI tools on synthetic plotting fixtures.

A passing setup probe establishes rendering capability at check time. It
does not establish that a later sample produced every plot: NanoPlot can
report a static-export problem while still producing HTML. Inspect the
specific run's outputs and logs. Full WSL installation, resource behavior,
scientific performance and representative-sample validation remain separate
acceptance checks.
