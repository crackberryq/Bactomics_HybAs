# AWARE+ — experimental development

The first implemented stage preserves taxonomic read groups from an existing
HybAs 9.0.2 identification-review receipt. The production HybAs entry point and
Bactologix interface are unchanged. This development branch is not a published
AWARE+ reconstruction release.

## Implemented: read preservation

`AwarePlus.smk` is a separate Snakemake entry point. Its required configuration
keys are `aware_plus_review` (the existing `review.json`) and `aware_plus_output`
(a new result directory). It reuses the controller environment and consumes
recorded merged reads and Kraken outputs; it performs no classification search.
The Python entry point is `workflow/scripts/aware_plus_preserve.py`, with
`--review` and `--output` arguments.

The receipt must declare ONT reads, and nonempty long-read data must actually be
present. Illumina-only inputs are rejected. PacBio is not supported by this first
entry point because the released HybAs receipt currently identifies ONT inputs.

All reads are retained in exactly one group:

| Group | Interpretation |
|---|---|
| `sample_species_name_taxid_N` | Kraken species assignment, including descendants, without detected conflicting species hits in the available report taxonomy |
| `_ambiguous` | Above-species assignments, conflicting species hits, missing hit taxonomy, or no informative hits |
| `_unclassified` | Kraken unclassified reads |

These are **provisional taxonomic groups**, not confirmed separate isolates.
Even a species-labelled group can contain related strains or classification errors.
The report explicitly marks linkage as unevaluated and reconstruction unavailable.
An individual species group with no assigned ONT reads is labelled accordingly.

Illumina mates remain paired. Out-of-order Kraken records are matched through a
disk-backed ID index. Duplicate IDs, unmatched reads, mismatched mates or lengths,
corrupt FASTQ, report/assignment count disagreements, changed evidence, and
conflicting taxonomic hierarchies stop the export. Existing result directories
are never reused by the exporter. A failed export removes only its own temporary
directory. A per-destination lock prevents concurrent writers; after a hard kill,
inspect the process and destination before removing a stale lock.

## Output bundle

- Each populated group contains `ont.fastq.gz` and/or paired
  `illumina_R1.fastq.gz` and `illumina_R2.fastq.gz`. Empty groups have no FASTQ files.
- `assignments.tsv.gz` records each read/pair, its group, the routing reason,
  candidate species associations and unresolved taxonomy. Candidate associations
  are not exhaustive when a taxon is absent from the available report hierarchy.
- `manifest.json` records input/output counts, checksums, the parent review digest
  and context, component source digest, runtime, elapsed/CPU time and peak memory.
- `report.html` displays the inventory and its interpretation limits.
- Snakemake retains the stage log and benchmark outside the protected bundle.

Shared/ambiguous reads are stored once. Their candidate associations allow later
analysis to consider them for multiple groups without treating those reads as
independent evidence or silently duplicating them into supposedly distinct genomes.

Only standard nucleotide Kraken output and six-column taxonomy reports are
accepted. Taxonomy absent from these reports remains unresolved. Validation is
currently software-only, using artificial reads and taxonomy labels. Read
partitioning accuracy, runtime and storage costs require mixed-sample benchmarking.

## Following stages and interface

1. Add the opt-in Bactologix AWARE+ panel to the identification review. It should
   show provisional groups, assigned long-read counts, unresolved reads and export
   status; short-read-only samples remain ineligible.
2. Develop and benchmark a long-read linkage assessment that records supporting,
   conflicting and non-informative evidence. Similarity alone is insufficient to
   assert that fragments came from the same organism.
3. Enable separately versioned reconstruction only after that assessment has a
   validated acceptance contract. Keep a supported reconstruction and a clearly
   labelled exploratory reconstruction as separate outputs. Unresolved sequence
   must remain traceable; neither output should be labelled guaranteed.

References: [Kraken 2 output formats](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats)
and [metaFlye strain complexity and long-read assembly](https://doi.org/10.1038/s41592-020-00971-x).
