# HybAs 10.0.0a9 development update

New runs with BUSCO and polishing stages generate a final-coordinate validation report and include its summary in the HybAs research report. For a completed run with an existing trajectory, launch this from Ubuntu inside a checkout of this revision:

```bash
bash scripts/refresh-polishing-validation.sh /path/to/run/outputs/isolate1 14
```

The first refresh installs a separate versioned analysis environment. It needs a network connection and free space for final-reference Illumina alignments. Intermediate alignment evidence is retained. The refresh inserts the summary into the existing research HTML if present. Its standalone companion is reports/polishing_validation/index.html. Failed attempts remain in validation-attempt-* directories. The existing trajectory is an input.

Stage-to-final comparisons use dnadiff one-to-one alignments. Adjacent-stage substituted and deleted bases use the earlier stage coordinates, while inserted bases use the later stage coordinates. Only bases uniquely projected onto the final sequence contribute to window totals. Affected bases across successive transitions can revisit the same locus. Windows are 10,000 bases per contig. Up to 20 nonempty windows are ranked by affected bases, breaking ties by contig and position. This is the rule for this release. Annotation overlap uses the union of final GFF rRNA intervals after exact sequence identity checking. It reflects the annotation source used in that run. Missing annotation or Illumina input is recorded explicitly.

Available trimmed paired Illumina reads are remapped against the actual final FASTA. Pileups report per-base evidence at projected positions and preceding anchors. Insertions and deletions are reported as anchored alleles. Mapping quality and base quality cutoffs are 20. These are observations from reads that can include data used for polishing. They are not an independent estimate of consensus accuracy.
