# HybAs AWARE 10.0.0a1

**Alpha target — unreleased.** Development branch:
`feature/v10-multiplatform-inputs`. Planned tag: `v10.0.0a1` (not published).
Work remains in [draft PR #4](https://github.com/crackberryq/Bactomics_HybAs/pull/4).

## Included development changes

- Resolve supported settings with Auto defaults and explicit overrides. Save
  the choices, reasons and sample plans before execution and reuse unchanged
  requests on resume. Detect local CPU/memory ceilings and expose the settings
  in computational records and reports. See [AUTO_SETTINGS.md](AUTO_SETTINGS.md).
- Pin the controller to Snakemake 9.27.1 and support its resource-value wrapper
  in the existing memory-budget check.
- Validate declarations for seven proposed input modes while preserving
  sample/library/run/read-group identity. This is metadata validation only;
  outputs explicitly retain `execution_ready: false`.
- Add an experimental FR paired-end library adapter: namespace/audit read
  pairs, preserve preparations through trimming, use separate SPAdes dataset
  entries, and record assembly/library provenance. SPAdes 4.3.0 has a separate
  pinned environment. Mate-pair and PacBio execution remain pending. See
  [ILLUMINA_LIBRARIES.md](ILLUMINA_LIBRARIES.md).
- Collect computational records automatically: all concrete rules have
  benchmarks, and each execution preserves job attempts, available CPU/memory/
  I/O measurements, timing, hardware and runner details, and JSON/TSV exports.
  Include computational history in Snakemake reports. Missing measurements
  remain unavailable; previous measurements are not counted as new work.

## Update and test

There is no published 10.0.0a1 update to install yet. Legacy execution routes
and runtime labels retain the v9 baseline. The experimental manifest adapter
is enabled explicitly through `input_manifest` in the run configuration; see
its guide for scope and limitations.

Before publication, align runtime/citation/version metadata, complete the
release checks, and register exact-commit Bactologix compatibility. The final
release notes must provide verified installation and migration instructions.
Existing runs retain their frozen configurations and source bindings.

## Verification and release scope

The foundation before Auto settings, at `27e9f9bd0f7ad17e60752f718b24a0b2fdb840a8`, passed
[197 tests on each of Python 3.11 and 3.12](https://github.com/crackberryq/Bactomics_HybAs/actions/runs/40417493724),
with native benchmarks enabled, plus
[controller and rule-environment deployment checks](https://github.com/crackberryq/Bactomics_HybAs/actions/runs/40417493819).
Auto settings adds capacity/override/snapshot tests and real Snakemake tests
with synthetic tool stubs for Auto execution, manual changes, resume, missing
databases, mode conflicts and setup previews. Current-head validation is
recorded on [PR #4](https://github.com/crackberryq/Bactomics_HybAs/pull/4).
These are software checks; they do not qualify new scientific routes.
This document is a development status record, not a release acceptance receipt.

Still pending:

- Scientific qualification of the FR paired-end adapter, plus executable
  adapters and qualification for the remaining proposed input modes.
- Integration of the separate [Research Report PR #3](https://github.com/crackberryq/Bactomics_HybAs/pull/3).
- Final Bactologix controls/compatibility, release metadata and publication.

See [V10_INPUTS.md](V10_INPUTS.md), [COMPUTATIONAL_DATA.md](COMPUTATIONAL_DATA.md)
and [RELEASE_CONVENTIONS.md](RELEASE_CONVENTIONS.md) for the detailed contracts.
