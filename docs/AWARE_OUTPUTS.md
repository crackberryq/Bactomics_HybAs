# Unified AWARE outputs (HybAs 10 development)

Regular AWARE remains the default. Target election, thresholds and per-platform cleaning decisions are unchanged. Extended preservation and statistical evidence are optional parts of the same module, not a separate analysis pipeline.

```yaml
kraken_decision_mode: aware
kraken_review_before_filtering: false
aware_run_until: assembly       # assembly | aware
aware_preserve_groups: false
aware_statistics: false
```

`aware_run_until: aware` completes identification, ordinary target cleaning and the requested AWARE outputs, then stops before trimming/assembly/polishing/annotation. `kraken_review_before_filtering: true` independently pauses earlier, after identification. Finalizing the review continues to the configured endpoint. Kraken must be enabled for AWARE outputs.

The output is `<sample>/reports/aware/`:

- `report.html`: primary and additional read-group inventory and optional statistics.
- `manifest.json`: `hybas-aware-bundle-v1`, relative file paths, SHA-256 checksums, supported assembly routes, parent context, original filtering history and continuation evidence.
- `primary/`: ordinary AWARE-cleaned inputs split back into their declared libraries/read groups.
- `groups/` when preservation is enabled: species assignment groups, one shared ambiguous group, one unclassified group and a read-assignment ledger. Each original read or pair appears exactly once in this preservation partition. The primary output is a separate view and can overlap these groups; they are not independent replicates.
- `filtering.tsv`: the unchanged native cleaning decisions.
- `statistics.json` when requested: original assignment-support evidence.

## Selecting and reusing outputs

Bactologix **AWARE outputs** lets a researcher select primary or provisional species groups and register them as derived samples. Each has its own normal run draft and configuration. **Samples → Import AWARE bundle** accepts the manifest after the entire output directory has been copied elsewhere. Original parent folders are not required to import it.

Each group has an embedded standard HybAs input manifest. Retained platforms, library preparation, pairing/orientation, read groups and PacBio HiFi/CLR chemistry determine its actual execution route. A hybrid parent can yield an Illumina-only child. Unsupported surviving combinations are listed with a reason rather than silently dropping a library. Legacy folder inputs are recorded as merged legacy libraries; declare an input manifest before the original run to retain more detailed library identities.

Ambiguous and unclassified groups are retained for inspection/export and are not automatically registered as isolates. Species assignments remain provisional: multiple strains, classification error and unresolved taxonomy may remain. This implementation does not evaluate long-read linkage or reconstruct strains. ONT is not required; Illumina, ONT and PacBio inputs are supported.

Derived Bactologix runs default to Kraken off so that the saved groups are not filtered again. The explicit **Reclassify derived inputs** setting applies the user's chosen Kraken settings again. Parent provenance remains in the sample/run record. Input bytes are checked on registration and before execution; lightweight metadata checks are used for UI readiness.

## Continue the parent assembly

Bactologix **Continue to assembly** creates a separate, checksummed `hybas-aware-continuation-v1` receipt, then validates the native continuation DAG and environment inventory. The original configuration stays frozen. Existing original/cleaned input evidence, context and bundle checksums must match. Explicitly requesting an assembly target cannot bypass a configured AWARE stop.

CLI users can create a JSON receipt containing `schema`, the absolute `bundle_path`, and its `bundle_sha256`, and supply its path and SHA-256 using `--config aware_assembly_continuation=... aware_assembly_continuation_sha256=...`. Review decision arguments, when applicable, must also be retained. Continuation of a stopped batch is deliberately one sample per run; create separate runs when interactive continuation is needed.

## Statistical interpretation

Statistics use original Kraken assignments, not selected/filtered groups. One Illumina pair is one sampling unit; each long read is one unit. The report gives 95% Wilson intervals for assignment support separately by platform and library/read group. They assume independent sampling units and are not calibrated probabilities of correct taxonomic identification.

Selection stability uses 200 deterministic multinomial bootstrap replicates (seed 1729), independently within each original library/read group. Counts are aggregated within each platform, regenerated as report counts and passed to the unchanged automatic target-election function. The report states the fraction of replicates selecting the original automatic recommendation, even when the researcher subsequently chooses a different manual target. It does not average platform support into a confidence score and does not change filtering decisions.

Database bias, duplicated molecules, correlated reads and sampling bias are not corrected by these statistics. Report count mismatches, duplicate identifiers and unresolved assignment taxonomy fail validation rather than produce unsupported statistics. Biological accuracy still requires external validation; software fixture tests do not establish it.
