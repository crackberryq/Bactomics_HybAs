# HybAs AWARE 10.0.0a6 — Alpha

This update is unreleased until the exact source passes all required checks and
is published as `v10.0.0a6`. See [Medaka model selection](MEDAKA_SELECTION.md)
for the pinned upstream evidence and complete policy.

## Changes

New configuration templates select `medaka_model_selection_policy:
documented_default`. A researcher-supplied model has priority, followed by an
exact mapping from one basecaller recorded on every read. When metadata is wholly
absent or identifies one unsupported basecaller, HybAs can use the documented
Medaka 2.2.0 consensus default `r1041_e82_400bps_sup_v5.2.0`. The installed Medaka
version, default identity and model-options source checksum must match the
reviewed upstream release.

The receipt and reports identify this as a fallback, with its reason, upstream
source commit, selected/resolved model and model SHA256. It is not a metadata
match or evidence of scientific suitability. An unmatched explicit model still
requires the recorded researcher acknowledgement, reason and exact model binding.

Mixed or partially annotated inputs require an explicit reviewed choice. Model
resolution, download, invalid-input, device and inference errors remain errors:
there is no retry using another model and no silent omission of Medaka. All-read
preflight, filtered-read reinspection and retained attempt outputs remain active.
Execution must agree with the recorded preflight policy and decision.

## Update and test

Use the new source and a new or revised draft. `strict` remains available, and an
existing configuration that omits the policy remains strict. Existing source
copies, saved settings, outputs and receipts are not rewritten. The new native
capability is version 2; Bactologix must approve the exact new source before
exposing the policy. Older approvals retain their original behavior.

Medaka stays at 2.2.0 and its environment is unchanged from 10.0.0a5. No assembly,
Racon, Polypolish, read-filtering or taxonomic thresholds are changed. Choosing a
fallback where a prior run stopped is a processing change, not a promise of exact
results from an earlier version or historical dataset.

## Verification and scope

Focused software tests cover precedence, strict compatibility, mixed/partial
metadata, runtime/source identity, errors without substitution, retained outputs,
preflight binding and report/Methods disclosure. The required native Medaka CI
check resolves absent, unsupported, matched and explicitly chosen cases against
bundled model weights and checks an invalid explicit choice fails without retry.
It downloads no additional model weights and runs no model inference.

Full regressions and the existing native environment, retention and rendering
checks remain publication requirements. Record the tested source commit and CI
URL in the published release. These checks establish software behavior, not
scientific accuracy, GPU inference behavior or equivalence to historical results.
