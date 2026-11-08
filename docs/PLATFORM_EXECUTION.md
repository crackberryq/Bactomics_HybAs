# HybAs AWARE 10.0.0a1 — experimental platform execution

This development alpha uses an explicit `input_manifest`; legacy folder inputs
retain Unicycler. Metadata validation is separate from software execution and
scientific qualification. A validator result of `execution_ready: false` is
intentional and is not a claim that the implemented adapter is unavailable.

| Declared inputs | Assembly | Optional polishing |
| --- | --- | --- |
| Illumina PE, repeated runs or complementary libraries | SPAdes 4.3.0, separate library entries | Polypolish |
| Illumina PE + supported mate-pair library | SPAdes, declared orientations | Polypolish from original PE libraries |
| ONT only | Unicycler | Racon → Medaka |
| Illumina + ONT | SPAdes hybrid dataset | Racon → Medaka → Polypolish |
| PacBio HiFi | Flye 2.9.6 `--pacbio-hifi` | Flye internal polishing |
| PacBio CLR | Flye 2.9.6 `--pacbio-raw` | Flye internal polishing |
| Illumina + PacBio HiFi or CLR | Flye long-read assembly | Polypolish from Illumina |

Auto selects the compatible stages; `none` disables the optional external
polishers. PacBio always uses Flye's one internal polishing iteration and never
receives ONT length/quality filters, Racon or Medaka. Explicit Racon/Medaka on
PacBio is rejected. CLR must be subreads with adapters/scraps removed and
individual passes separated, not unprocessed polymerase reads. All FASTQ must
be four-line records with Phred+33 quality encoding. Chemistry is declared,
never guessed from filenames. BAM input is not implemented.

## Illumina preparation

Each library declares `library_type` (`paired_end` or `mate_pair`), orientation
(`fr`, `rf`, `ff`) and read groups. Ordinary FR paired-end libraries default to
fastp. Other orientations and mate-pairs require `preprocessing`:

- `pretrimmed`: user declares adapters/junctions already removed. Preserve the
  declared orientation/type; do not run fastp or label a report as a fastp run.
- `nextera_mate_pair`: raw Nextera mate-pair RF input only. NxTrim 0.4.3 with
  `--rf` produces RF mate-pairs, FR paired-end, unknown, and singleton reads.
  Preserve the first two categories as distinct SPAdes datasets. Retain unknown
  and singleton reads as unpaired evidence. No additional fastp step is applied.

Declared insert sizes order libraries; SPAdes estimates insert distributions.
Mate-pair-derived reads are excluded from the pooled polishing inputs. Multiple
independent PE preparations require a declared repeated-run or complementary
strategy. Complementary libraries require insert sizes. Mixed ONT/PacBio,
mixed HiFi/CLR, mate-pair-only, and mate-pair + long-read inputs are rejected.
Specialized protocols beyond these declarations are not auto-detected.

## Audit and reporting

The workflow preserves original headers, namespaces IDs before AWARE, audits all
FASTQ records, and keeps library/read-group identities. Empty or inconsistent
inputs fail rather than silently dropping a library. AWARE target selection,
per-platform cleaning and confidence policies are unchanged. PacBio uses its
own classification, cleaning, QC and report labels.

`reports/hybas_research/hybas_report.html` is created by the workflow, with
SQLite-backed paired/long-read reconciliation, taxonomic removal tables, assembly
stages, and input plans. Audit assignment counts are classification evidence,
not proof of biological origin. When Kraken is off, assignment verification is
unavailable; read-retention counts still describe the files. Computational
history is recorded for every invocation and linked from the report. All new
concrete rules have benchmarks. Flye's memory resource is a scheduler reservation,
not an enforced process memory limit.

Mocked workflows check routing, failure/resume, no-op reuse and input provenance.
Native CI checks Flye chemistry flags and NxTrim categories on deterministic
random reads and SPAdes on a random two-library dataset. These are software
contracts; representative reference-data qualification remains required before
claiming analytical performance.

Primary contracts: [Flye 2.9.6 usage](https://github.com/mikolmogorov/Flye/blob/2.9.6/docs/USAGE.md),
[NxTrim](https://github.com/sequencing/NxTrim), and
[SPAdes datasets](https://ablab.github.io/spades/running.html#specifying-multiple-libraries-with-yaml-data-set-file).
