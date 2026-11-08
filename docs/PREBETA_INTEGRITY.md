# Pre-beta evidence integrity corrections

This development change closes the computational/reporting defects identified in
H1–H5 of the 8 November 2026 review. It does not qualify experimental routes or
establish assembly accuracy, taxonomic accuracy, or environmental suitability.

## Medaka execution and model provenance

Every scheduled Medaka execution uses a new retained attempt directory. Native
alignment/inference caches from an earlier execution cannot satisfy the stage.
Normal unchanged Snakemake reruns remain no-ops; a forced execution recomputes.
The `current/consensus.fasta` publication path leaves a pre-upgrade
`medaka/consensus.fasta` intact. New completed and failed attempt records are
retained under `current/attempts/`; these consume storage and may be removed only
as an explicit workspace maintenance decision after preserving needed evidence.

ONT read namespacing preserves comments, including basecaller model tags and
Dorado read groups. Automatic selection requires one resolved basecaller identity
on every input read, checks the whole FASTQ rather than just the first 100 reads,
and invokes pinned Medaka's model resolver without a default fallback. Missing,
incomplete or mixed metadata requires researcher review and explicit
`medaka_model` configuration. This changes the previous permissive behavior:
an unlabelled FASTQ no longer silently receives Medaka's default model.
An explicit model still passes Medaka's native compatibility check; its selection
is recorded as explicit, not metadata-verified. Do not infer model suitability
solely because the tool completed.

The receipt records input hashes, the effective resolved model file/hash, metadata
coverage, command, native attempt directory and consensus hash. Inputs/model are
checked again before publication. Consensus and receipt each use atomic file
replacement; they are not a two-file transaction. Downstream promotion checks
that receipt status and checksum match the consensus. Configured model file
content changes invalidate the stage. Named/default model stores are resolved
and hashed when execution is scheduled; provenance does not promise continuous
monitoring of mutable installed environments during a no-op rerun.

## Processing checks and scientific interpretation

The former `overall_qc_status` output is replaced by
`processing_evidence_status` (`RECORDED`, `REVIEW_REQUIRED`, `MISSING`) and a
separate `analytical_suitability=NOT_ASSESSED`. BUSCO 0% can be a valid recorded
measurement and is never converted into an analytical pass. Values outside
0–100% or nonfinite values receive a recording warning. No new scientific
acceptance threshold has been introduced. HTML makes the interpretation visible
in the overview and quality section and explains legacy bookkeeping PASS values.
MultiQC custom-content headings use the new separate status fields.

## Taxonomic confidence

The legacy confidence path now compares alternative taxa at the selected rank,
so a selected species' own genus is not treated as a competing population.
Target selection and independent per-platform cleaning are unchanged. Existing
single-platform handling remains unchanged. Confidence is still a heuristic,
not a calibrated probability of identity.

## BUSCO resource identity

The prepared lineage's file content hashes (excluding the readiness marker)
produce a deterministic resource digest. The digest is a scheduler parameter and
all resource files are declared inputs. The lineage is checked before and after
BUSCO execution and a receipt is saved. Modified contents invalidate BUSCO even
if path, file size, timestamp and `.complete` are unchanged. Content receipts and
file hashes are retained in provenance and the summary records the executed
dataset digest. Datasets must remain immutable while analysis is running.

## Verification scope

Regression tests use harmless synthetic FASTQ/FASTA, fabricated model files and
native-command stubs to verify identity, publication, rerun and presentation
behavior. They do not test neural inference or BUSCO biological accuracy.
Real pinned-tool acceptance, representative reference-data benchmarking and
clean-install qualification remain required before a supported research beta.

The command contract was checked against the pinned upstream
[Medaka 2.2.0 consensus wrapper](https://github.com/nanoporetech/medaka/blob/v2.2.0/scripts/medaka_consensus)
and [model metadata parser](https://github.com/nanoporetech/medaka/blob/v2.2.0/medaka/models.py).
Local checks exercised the real Snakemake scheduler with tool-contract stubs;
the documented `HYBAS_TEST_DISABLE_BENCHMARKS=1` hosted-environment fallback
excludes resource sampling only in temporary test copies. Production benchmark
code remains unchanged and is not validated by that fallback.
