> This page describes the original FR adapter. See [PLATFORM_EXECUTION.md](PLATFORM_EXECUTION.md) for the expanded execution matrix and preparation options.

# HybAs AWARE 10.0.0a1 — Illumina library execution

**Experimental development adapter, unreleased.** This is software integration,
not qualification of assembly accuracy. It belongs to draft PR #4 and the
[10.0.0a1 alpha target](RELEASE_10.0.0a1.md).

## Executable scope

| Declared input | Execution in this adapter |
|---|---|
| One FR paired-end preparation | SPAdes with one library |
| Repeated lanes/runs of one FR paired-end preparation | Preserve read groups, trim together within that library |
| Independent FR paired-end preparations | Trim separately; separate SPAdes dataset entries |
| Complementary FR paired-end preparations with distinct declared insert sizes | Separate dataset entries ordered by declared insert size |
| Mate-pair, RF or FF orientation | Explicit preprocessing required; see expanded guide |
| ONT or PacBio manifest modes | Supported by separate platform adapters; see expanded guide |

Legacy folder/sample-sheet runs continue to use the existing Unicycler path.
Do not put independent preparations into a legacy raw folder and assume their
library boundaries will be inferred. This adapter requires explicit metadata.
All libraries must be ordinary paired-end genomic DNA inputs using Phred+33
quality scores. Legacy Phred+64 input needs conversion outside this adapter;
quality encoding is not inferred from filenames or uniformly high scores. A preparation kit
name does not automatically select MDA, RNA, UMI, linked-read or other special
processing; those protocols are outside this adapter's contract.

## Configuration

Start from `config.example.yaml` in a new development run directory. Set:

```yaml
input_manifest: /absolute/path/to/inputs.local.yaml
sample_sheet: ''
run_mode: auto
```

The manifest supplies sample IDs. An optional sample sheet can still provide
taxon, lineage and model overrides, but must contain exactly the same sample IDs.
Use `examples/v10/illumina-libraries.example.yaml` as the declaration template;
its paths are fictional. Read paths resolve relative to the manifest directory.
The manifest itself resolves from the working directory. Explicitly choose
`dual_illumina_runs` or `dual_illumina_complementary` for independent preparations.
Auto refuses to guess their intended relationship. The existing top-level
`run_mode` remains `auto` or `illumina_only` for these inputs; the more specific
strategy is saved separately in the library plan and run summary.

Normal Snakemake execution with `--use-conda` installs the pinned
`library_assembly` environment (SPAdes 4.3.0) when selected. Existing trimming,
QC and optional polishing environments retain their own pins. No direct
assembler override or arbitrary extra flags are introduced.

## Data and decision contracts

1. Validate the whole manifest and executable scope before starting jobs.
   Check file existence/initial FASTQ structure during preflight; audit every
   paired FASTQ record during preparation before classification. Fail on unequal
   counts, mismatched mate IDs, wrong mate labels or malformed records.
2. Namespace read IDs by library, read group and pair number. Preserve original
   headers in `work/inputs_merged/library_identity/original_headers.jsonl.gz`,
   with the declared plan and read-group pair counts alongside it. Raw inputs
   are unchanged. Namespacing prevents repeated sequencer IDs from crossing
   library boundaries during extraction.
3. Classify the combined Illumina stream once and apply the existing AWARE
   sample/platform decision unchanged. Library metadata does not change target
   election, confidence or cleaning policy.
4. Restore library boundaries from the retained IDs, then run fastp separately
   per preparation. Refuse unknown/reordered identities or a completely empty
   library after cleaning/trimming; never silently omit a declared preparation.
5. Pass each trimmed library as a separate paired-end/FR SPAdes dataset entry.
   Declared insert sizes are retained and used for ordering; SPAdes estimates
   insert distributions. No unsupported insert-size dataset field is invented.
   SPAdes receives explicit `--phred-offset 33`, matching fastp's output encoding.
   This is standard SPAdes execution, without automatic special-purpose modes.
6. Promote `contigs.fasta` through `min_fasta_length` to `stage_spades.fasta`.
   Existing optional short-read polishing and downstream QC/reporting consume
   the stage. Combined trimmed streams remain available for aggregate QC and
   independent-mate Polypolish mapping; assembly receives the separate libraries.

`trimmed/libraries/` contains trimmed pairs, individual fastp reports,
`dataset.yaml` and `libraries.json` recording dataset order and counts. The
sample fastp HTML embeds the individual reports. Its JSON uses
`hybas-library-fastp-v1`, not the single-library fastp schema; it preserves each
report without inventing pooled percentages. MultiQC can discover the individual
reports. Consumers parsing the old fastp JSON must use the per-library files.

Run summaries identify `assembler: spades`, the strategy and library count.
Polishing provenance records the real SPAdes source stage, including when
polishing is disabled or skipped. The separate environment export and SPAdes
version are captured. The frozen settings identity includes the manifest hash;
its per-sample plan includes declared libraries and selection. A changed
manifest creates a new settings record. Unchanged requests reuse their record.
Changing only polishing settings reuses completed preparation, classification,
trimming and assembly. All adapter jobs have computational benchmarks and appear
in the existing per-invocation records, including failures and retries.

## Validation and remaining work

Synthetic tests exercise library/read-group separation, duplicate original IDs,
spaces in paths, complete pair audits, post-cleaning separation, blocked modes,
minimum-length promotion, optional polishing, failure/resume and no-op reuse.
CI also installs the pinned environment and runs a two-library SPAdes smoke test
on deterministic random toy reads. These checks establish software contracts,
not biological accuracy or suitability for a particular organism/preparation.

Representative benign reference-data qualification remains a release gate.
Mate-pair junction processing, other declared orientations and PacBio HiFi/CLR
adapters are now implemented; see the expanded guide. Special preparation
protocols and release registration require separate qualification. The standalone manifest validator deliberately still reports
`execution_ready: false`; metadata validity alone never certifies an executable
or scientifically qualified route.

Primary tool contracts: [SPAdes 4.3.0](https://github.com/ablab/spades/releases/tag/v4.3.0),
[SPAdes dataset format](https://ablab.github.io/spades/running.html#specifying-multiple-libraries-with-yaml-data-set-file).
