## 10.0.0a4 — explicit swap scheduling

- Honor `bactologix_compute.swap_mb` as an additional scheduler allowance when the supplied Snakemake `mem_mb` contains both RAM and swap. Keep native physical-RAM ceilings, the 80% Auto policy, and all tool memory defaults separate and unchanged.
- Record the requested RAM/swap split and effective scheduler budget in resolved-settings provenance; show the effective split in the controller log. Existing resolved records are retained unchanged.
- Reject malformed allowances and a combined budget without a positive RAM portion. Retain the database-size reservation guard and exercise it with pinned Snakemake 9.27.1 dry runs.
- Usable swap must be verified by the launcher. Bactologix performs current free-swap and container-limit checks; the standalone workflow does not measure or reserve swap.

## 10.0.0a3 — native output retention

Preserve BUSCO 6.0.0 Python-managed cleanup products and separate stdout/stderr for every native job, including failed jobs. Preserve Flye 2.9.6 consensus/polishing intermediates and previous logs/parameter snapshots before its reviewed Python cleanup. Source-checked launchers retain normal cleanup semantics and analysis arguments; retention failures stop the run. Native package contract tests gate publication. See [output retention](docs/OUTPUT_RETENTION.md) for exact scope and exclusions.

## 10.0.0a2 — Conda controller correction

Use the executing Snakemake module version for retention compatibility and run version records. Accept known Conda distribution placeholders only with a supported runtime; reject contradictory versions. Native Conda controller checks now gate release publication. No analysis parameters change.

## 10.0.0a2 — output retention

Retain prior job products, native tool sidecars, logs and benchmarks before Snakemake removes outputs for reruns or retries. Failed and temporary products are retained. Enable supported native keep options and retain generated SAM/FASTQ streams. No polishing, filtering or assembly-selection parameters change. See [output retention](docs/OUTPUT_RETENTION.md) for paths, storage implications and upstream limits.

## HybAs AWARE 10.0.0a1 — integration update (unreleased)

- Integrate regular AWARE 9.0.2 identification review and the branded report into the v10 development branch.
- Pause before cleaning/assembly for legacy and manifest inputs, including PacBio; preserve platform labels in evidence and reports.
- Reuse frozen Auto resource settings on continuation while recording the reviewed decision separately. Include manifest content and raw PacBio identity in stale-review checks.
- Keep AWARE+ separate. No production release or compatibility registry activation.

## HybAs AWARE 9.0.1

## HybAs AWARE 9.0.2 — Alpha

- Optional identification review stops after Kraken and emits an early HTML/JSON report.
- Manual targets can be deferred; a recorded, checksummed decision resumes the existing work.
- Changed evidence, inputs, database identity or workflow settings reject stale continuation.
- The existing target election, platform cleaning and confidence rules are unchanged.

The primary assembly report now uses the Bactomics visual design: navy sidebar, assembly overview, separate read-retention panels, AWARE target and platform-cleaning evidence, recorded polishing route, genome quality, annotation counts and methods/provenance.

New runs write `reports/hybas/report.html`, `reports/hybas/results.json` (schema `hybas-report-1`) and `reports/hybas/methods.txt` per isolate. The report rule depends on existing analysis outputs and executes no analysis tools. MultiQC remains available in its original location as a detailed tool report.

The standalone renderer accepts the saved bundle or a legacy `metadata/run_summary.json`. Bactologix a81 exposes report-only regeneration using this installed renderer; it saves a separate edition and preserves original results. Legacy summaries display available metrics and mark uncaptured details unavailable.

No AWARE selection, cleaning or confidence policy is changed. No circularity or stage-quality improvement is inferred where it was not recorded. Read-count charts keep Illumina mates and ONT separate. BUSCO uses the recorded lineage; coverage remains explicitly a base-count estimate. Downloads include saved report evidence and a methods draft for researcher review.

# Changelog

## HybAs AWARE 10.0.0a1 — unreleased alpha target

- Add an experimental FR paired-end library execution adapter. Audit and
  namespace read pairs, keep preparations separate through fastp and SPAdes
  4.3.0, and retain library/stage provenance and computational benchmarks.
  Manifest content participates in frozen settings identity.
- Extend execution to explicit pretrimmed RF/FF libraries, Nextera mate-pairs,
  ONT manifests and separate PacBio HiFi/CLR paths. Preserve library identity
  and platform-specific evidence; scientific qualification remains pending.
- Generate the integrated research report with read-retention audits, actual
  platform labels, input plans and computational-history links. See
  `docs/PLATFORM_EXECUTION.md` for the supported input/preparation matrix.

- Resolve supported run controls with Auto defaults while retaining explicit
  overrides. Detect local CPU/memory limits, verify manual input modes and
  retain existing processing presets and required-database checks.
- Freeze resolved choices/reasons before execution, reuse them for unchanged
  requests and retain earlier snapshots after explicit configuration edits.
  Include selections in computational records, reports and provenance.
- Adopt BaSIS-style product names, numeric prerelease versions, `v<version>`
  tags and version-specific release documents for future updates.
- Record completed work, pending execution/release integration and
  release gates in `docs/RELEASE_10.0.0a1.md`. This documentation update does
  not publish a release or change the v9 runtime/configuration contract.

## v10 development — computational records

- Enable extended native benchmarks for every concrete rule, including QC and provenance stages that previously had no benchmarks.
- Persist a separate computational record for each execution invocation, including hardware/runner details, requested resources, UTC timestamps, job attempts, available resource measurements and raw benchmark archives.
- Preserve failures, retries, resumes and no-work invocations; never reuse an old benchmark as a new measurement. Missing metrics remain null/NA.
- Include the computational history in Snakemake reports and provide JSON/TSV exports. Run elapsed time is separate from summed CPU time and the largest individual job RSS peak.
- Bactologix 0.1.0a17 automatic release-specific runner preparation is now published.

## v10 development — runner compatibility

- Pin the development runner and regression CI to Snakemake 9.27.1, retaining Python 3.11 in the Conda recipe.
- Accept newer Snakemake resource objects in the existing memory-budget check, retaining the earlier scalar representation and the same reservation limits.
- Document release-specific automatic runner preparation in Bactologix 0.1.0a17. New multi-platform manifests remain metadata-only until their execution adapters are qualified.

## Unreleased — package layout only

- Group rule environments, report captions, helper code, scripts and MultiQC configuration under `workflow/`.
- Group the sample template and CI example inputs under `examples/`.
- Keep one canonical citation file and move the license and contribution/community files to standard repository locations.
- Update file references, CI paths and documentation. The workflow version, analysis settings, stages, commands and isolate output layout are unchanged.

## 9.0.0-AWARE — 2026-10-19

Engineering update from v8.8 commit `b14ff7688b72f1ae4a39bd99a4305d16a0e53ec1`.

### Correctness

- Enable Illumina-only Polypolish under the existing `tripolish` interface.
- Promote Medaka when hybrid Tripolish disables Polypolish.
- Supply two independent, all-alignment SAM files to Polypolish; retain secondary records with `SEQ=*`. Preserve separate paired mapping outputs.
- Record requested/planned/executed/skipped stages, final source and FASTA hashes; validate stage FASTA before promotion.
- Reject malformed manual taxids, duplicate/unsafe isolate names and unmatched lane pairs before processing. YAML null remains empty.
- Remove shadowed AWARE helpers while preserving the formerly active default behavior. Honor all four configured confidence thresholds.
- Parse native QUAST TSV orientation. Make MultiQC depend on the generated summary and configuration. Record disabled QC modules as skipped.

### Execution and provenance

- Declare assembly, polishing and database-aware Kraken memory reservations; serialize Kraken by default. Forward a configurable SPAdes memory setting and minimum export length.
- Bound compression threads and separate BWA from sorting to avoid simultaneous full-core allocation.
- Preserve assembly diagnostics and intermediates on failure; remove polishing's mutation of an input BAM.
- Track execution-affecting configuration in rule parameters without making reporting changes invalidate classifications.
- Capture source version/hash separately from the configured version label; track effective configuration and optional environment manifests.
- Validate BUSCO dataset structure, use local datasets offline during normal runs, and preserve incomplete data during explicit preparation.
- Restore report captions and Conda configuration; provide a lean native Snakemake controller environment, deployment checks and portable examples.
- Retain normal Snakemake stage messages and write detailed tool output to durable logs.

### Compatibility and validation limits

Default AWARE target election, independent platform actions and conservative confidence are protected by frozen-function regression tests. Nondefault confidence thresholds now affect labels as documented.

Polished FASTA may change because the Polypolish input defect is corrected. Illumina-only output can now be polished. Do not assume byte-for-byte assembly equivalence with v8.8. Keep old checkouts, configurations and outputs for comparison.

The new software tests use synthetic fixtures and mocked external tool contracts. They test orchestration, not assembly accuracy. Historical phase 1–3 validation assets are not bundled or recreated. Real-tool biological comparisons, persistent Kraken reuse and new assembly-selection policies are outside the completed validation scope.

