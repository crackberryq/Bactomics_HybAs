# HybAs 10 development: input metadata foundation

Status: **development preview, manifest validation only**. The scientific
routes remain the HybAs 9.0.0-AWARE baseline, with Snakemake 9.27.1 pinned for
this development branch. This change does not implement a new assembler,
change any AWARE decision, or connect these manifests to `Hybas.smk`.
The Research Report continues independently in PR #3.

## What works now

The standalone validator accepts a versioned JSON or YAML manifest, checks
metadata consistency across the entire batch, and emits normalized JSON to
standard output. It retains sample, library, sequencing-run and read-group
identities. It accepts declarations for all seven proposed modes, with separate
PacBio HiFi/CLR metadata. Every result includes `execution_ready: false`.

```bash
python workflow/scripts/validate_input_manifest.py examples/v10/inputs.example.yaml
```

The example contains fictional paths. The validator does not open sequence
files, check whether they exist, audit read pairs, or assess their quality. A
successful exit validates **declarations**, not the data or assembly route.
It reads only the manifest and writes only to standard output. Exit code 0
means valid metadata; code 2 means a manifest/read error, with the field location
on standard error. No Snakemake or external analysis tools are launched.

Copy the example to `inputs.local.yaml` before editing. Local YAML files are
already gitignored. Do not replace the existing v9 `config.yaml` or
`samples.tsv` with this manifest; v9 does not consume it.

## Schema version 1

Unknown keys are errors at every level. Duplicate JSON/YAML keys are errors.
Identifiers start with a letter or digit and then use letters, digits, dots,
underscores or hyphens. Sample IDs are unique in a manifest; library IDs are
unique within a sample; read-group IDs are unique within a library.

| Object | Required fields | Optional fields |
|---|---|---|
| Manifest | `schema_version: 1`, nonempty `samples` list | None |
| Sample | `sample_id`, nonempty `libraries` list | `mode` (default `auto`) |
| Illumina library | `library_id`, `platform: illumina`, `library_type`, `orientation`, `read_groups` | `insert_size_bp` |
| ONT library | `library_id`, `platform: ont`, `read_groups` | None |
| PacBio library | `library_id`, `platform: pacbio`, `chemistry`, `read_groups` | None |
| Illumina read group | `read_group_id`, `run_id`, `r1`, `r2` | `lane_id` |
| Long-read group | `read_group_id`, `run_id`, `reads` | `lane_id` |

- `library_type`: `paired_end` or `mate_pair`.
- `orientation`: explicitly declare `fr`, `rf`, or `ff`; never inferred from filenames.
- `insert_size_bp`: positive integer representing the declared library insert
  size. Required for each complementary Illumina library; it is not measured by
  this validator or forwarded to an assembler.
- PacBio `chemistry`: explicitly declare `hifi` or `clr`. Corrected CLR, BAM
  imports and other representations need a later schema/adapter decision.
- All files must be literal Linux/WSL paths ending in `.fastq`, `.fq`,
  `.fastq.gz`, or `.fq.gz`. Paths are resolved relative to the manifest's
  directory. No globbing, downloading, Windows-path translation or read
  conversion is performed. Spaces inside paths are allowed.
- A file may appear only once across the whole batch. Normalized relative,
  absolute and existing symlink aliases are checked; content identity and
  hard-link aliases are not audited.

Put multiple lanes/runs of the **same physical library** under that library's
`read_groups`. Declare independent library preparations as separate libraries,
even when the files came from the same sequencing run. This validator never
merges or flattens library groups.

## Mode selection contract

Here `auto` means **infer a metadata mode from explicit declarations**. It does
not automatically select an assembler, polishing algorithm or executable route.
Later execution integration must enforce qualification separately.

| Mode | Accepted declarations |
|---|---|
| `illumina_only` | One paired-end library, one sequencing run; multiple lanes allowed |
| `dual_illumina_runs` | One paired-end library with at least two run IDs, or at least two independent paired-end libraries |
| `dual_illumina_complementary` | Paired-end plus mate-pair libraries, or at least two paired-end libraries with distinct declared insert sizes |
| `ont_only` | ONT libraries only |
| `illumina_ont` | Paired-end Illumina library/libraries plus ONT |
| `pacbio_only` | PacBio libraries of one explicitly declared chemistry |
| `illumina_pacbio` | Paired-end Illumina library/libraries plus one PacBio chemistry |

The names `dual_illumina_*` describe strategies; the schema allows more than two
read groups/libraries. With two or more paired-end libraries and no long reads,
`auto` refuses to guess between extra coverage and complementary-library intent.
Choose the mode explicitly. Complementary declarations require insert sizes;
all-paired-end complementary libraries must have at least two distinct sizes.

Mixed ONT/PacBio, mixed HiFi/CLR, mate-pair-only, and mate-pair plus long-read
combinations are not supported in this first contract. A manual mode cannot
silently omit incompatible input. Library orientation and run IDs are metadata
claims; checking them against the actual files is future work.

The output preserves the original requested mode, the selected metadata mode,
the reason category (`explicit` or `declared_metadata`), each library/read group,
and a SHA-256 of the source manifest. This hash identifies the manifest bytes,
not the sequence files. Input declaration order is retained.

## Assembly/polishing design matrix — not executable

These are candidate integration directions, not tested v10 tool selections.
No new scientific tool versions or assembler commands are introduced in this milestone.

| Mode | Assembly direction | Polishing direction / open decision |
|---|---|---|
| Illumina only | Preserve existing Unicycler route | Preserve configurable short-read polishing |
| Dual Illumina, repeat runs | Preserve same-library grouping; keep independent libraries separate; evaluate multi-library SPAdes adapter | Library-aware short-read evidence and provenance |
| Dual Illumina, complementary | Evaluate multi-library SPAdes adapter | Preserve library type/orientation through preprocessing and polishing eligibility |
| ONT only | Preserve existing Unicycler baseline; evaluate a separate long-read assembler adapter | Explicit ONT metadata/model compatibility before any new route |
| Illumina + ONT | Preserve current hybrid baseline | Preserve existing configurable polishing contract |
| PacBio only | Separate HiFi/CLR adapters; Flye is a candidate | Chemistry-specific qualification; no inherited ONT polishing defaults |
| PacBio + Illumina | Qualify HiFi and CLR combinations independently | Establish whether/how short-read polishing is appropriate for each route |

Primary references: [SPAdes input model](https://ablab.github.io/spades/input.html),
[SPAdes dataset representation](https://ablab.github.io/spades/running.html#specifying-multiple-libraries-with-yaml-data-set-file),
[Flye supported input](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md),
and [Unicycler](https://github.com/rrwick/Unicycler).
Upstream support does not establish that HybAs's installed versions, wrappers or
postprocessing support a combination. This matrix intentionally does not claim
that any candidate has been qualified here.

## Development gates

1. **This milestone:** input declarations, strict validation, metadata-only mode
   selection, examples, and regression tests.
2. **Execution adapters:** define stable library-level artifacts, provenance and
   resource contracts; add independent adapter tests before DAG integration.
3. **AWARE extension:** keep target election, independent actions and confidence
   separate. New platforms must have explicit evidence identities and a
   separately reviewed policy; PacBio must not masquerade as ONT.
4. **Qualification:** real-tool compatibility and representative benign
   reference-data benchmarks for each proposed route; compare with the frozen
   baseline. Software tests alone do not establish accuracy.
5. **Release integration:** report/UI controls, capability-gated execution,
   migration documentation and a release from the tested commit.

## Runner installation and compatibility

This branch's `environment.yaml` pins Python 3.11 and Snakemake 9.27.1. CI uses
that same Snakemake release on Python 3.11/3.12. A two-line compatibility change
unwraps Snakemake's newer `Resource.value` object in the memory-budget check;
it also accepts the scalar resource values used by 9.6.2. Reservation values,
AWARE decisions, scientific routes and rule environment recipes are unchanged.

The full local software suite passes 184 tests on 9.27.1 using the documented
container benchmark workaround; 13 existing DAG/budget checks also pass on
9.6.2. Native Conda installation and normal benchmarking remain CI checks.

Bactologix 0.1.0a17's proposed installer reads the selected release's environment
file, installs/checks that release's runner, then activates it. Different
requirements use a separate environment. Older runners and cached rule
environments are retained. No global "install latest Snakemake" operation is
used. The unfinished v10 source remains outside Bactologix's reviewed execution
catalogue until adapter integration and qualification are complete.

The production v9 branch remains unchanged. Metadata fixtures contain no
sequence data, and mocked workflow tests do not establish biological accuracy.
