# Changelog

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
