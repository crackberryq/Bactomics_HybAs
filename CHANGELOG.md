# Changelog

## HybAs AWARE 10.0.0a1 — unreleased alpha target

- Add an experimental FR paired-end library execution adapter. Audit and
  namespace read pairs, keep preparations separate through fastp and SPAdes
  4.3.0, and retain library/stage provenance and computational benchmarks.
  Manifest content participates in frozen settings identity. Mate-pair, RF/FF
  and PacBio execution remain blocked pending separate adapters/qualification.

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
