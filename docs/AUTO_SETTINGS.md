# Automatic run settings

HybAs AWARE 10.0.0a1 development resolves supported controls automatically by
default. An explicit YAML or Snakemake `--config` value overrides Auto. The
requested YAML is never rewritten. This applies to the existing Illumina-only,
ONT-only and hybrid routes; new manifest modes still require execution adapters.

## Controls

| Setting | Auto behavior | Manual override |
|---|---|---|
| `run_mode` | Detect each sample's existing input layout | `illumina_only`, `ont_only`, `hybrid`; incompatible input is an error |
| `threads` | Use CPU capacity within Snakemake's core ceiling | Positive integer; Snakemake can scale job threads to its ceiling |
| `total_mem_mb` | Use 80% of detected physical/container memory, bounded by a supplied Snakemake memory ceiling | Positive MB budget within the machine limit |
| `assembly_mem_mb` | Up to the existing 8192 MB reservation, bounded by the budget; enough for an explicit SPAdes limit | Positive MB reservation |
| `spades_memory_gb` | Fit within the assembly reservation using the existing 1024 MB/GB check | Positive integer |
| `polishing_mem_mb`, `busco_mem_mb`, `prokka_mem_mb` | Up to the existing 4096 MB reservations, bounded by the budget | Positive MB reservation |
| `kraken_mem_mb` | Existing database-size estimate plus overhead | Positive MB reservation |
| `polish_mode` | Existing Tripolish preset, with input-compatible stages selected per sample | `none`, `racon`, `medaka`, `tripolish` |
| `run_kraken`, `run_busco`, `run_prokka`, `run_multiqc`, `run_polypolish` | Existing enabled defaults and route eligibility; Kraken's explicit `off` policy disables its Auto flag | `true`/`on` or `false`/`off` |

Omitting a supported Auto control has the same meaning as `auto`. Explicit
numeric/boolean values in existing YAML remain manual values. Invalid values,
including null where a positive number is required, are errors, not Auto.
Conflicting explicit Kraken flags/policies retain the existing validation.

`kraken_decision_mode: auto` retains its original Kraken algorithm meaning.
It is not a request to choose between `aware`, `manual` and `off`. The template
continues to use `aware`. Scientific thresholds, retention percentages, Racon
rounds, database locations, lineage/model identifiers and manual taxids keep
their existing contracts. Auto does not infer them from filenames or alter
AWARE election, platform actions or confidence calculations.

Auto does not silently disable requested QC because a required database is
missing. Normal preflight still fails and explains what must be prepared.
An undersized machine likewise produces a capacity error when a required
reservation cannot fit; the Kraken estimate is never clipped to hide it.

## Capacity and scheduling

Detection targets the local Linux/WSL controller. It considers logical CPUs,
CPU affinity, cgroup v1/v2 CPU quotas/cpusets and physical/cgroup memory limits,
including ancestor limits. The default profile requests all available cores;
the workflow bounds the scheduler by detected capacity. Explicit `--cores`
and `--resources mem_mb=...` remain additional ceilings. If memory cannot be
detected, supply `total_mem_mb` or a Snakemake memory budget explicitly.

Memory budgets use decimal MB from detected bytes. Existing per-tool defaults
and the conservative SPAdes 1024 multiplier are retained. Budgeting is based
on capacity rather than momentary free memory, so unrelated running processes
do not change Auto choices on every launch. Leave additional headroom using a
manual budget when sharing the machine. These are scheduling reservations,
not operating-system limits or a guarantee of actual peak memory usage.
Remote-worker capacity and GPU/CUDA/model selection are outside this resolver.

## Preview, execution and resume

A normal dry run prints the selected CPU/memory limits and settings-record
path. It performs validation but creates no settings snapshot. Immediately
before execution, the start hook saves a complete, checksummed JSON record:

`<base_dir>/provenance/runtime/resolved_settings/<settings-id>.json`

The record contains the original request, resolved configuration, selection
reasons, capacity at resolution and each sample's detected mode, metadata,
input paths, enabled modules and planned polishing stages. Computational run
records reference this identity, keep their own JSON copy and show the choices
in a collapsible section of the computational report. Workflow provenance
also records the requested and effective configurations and snapshot path.

The same request, sample-sheet content, working directory and resolver policy
reuse the saved choices. Changes in machine capacity do not silently choose a
different preset. Current CPU/memory ceilings still apply, and a reservation
that no longer fits stops before processing. Changing only `--cores` can scale
job threads without replacing the saved requested thread count.

An explicit configuration edit or sample-sheet change creates another record;
earlier records remain untouched. To deliberately re-detect capacity for the
same YAML, set a new `settings_revision`, for example `settings_revision:
workstation-2`, or create a new draft. If input routes or input paths change
under an unchanged request, HybAs asks for a new revision rather than changing
the frozen plan. Missing required inputs still fail normal validation first.
If a database update changes the saved Auto Kraken memory estimate, a new
revision is also required so the old reservation is not silently reused.

Dry runs, touch operations and report-only rendering do not create snapshots.
`--no-hooks` is rejected for controller execution because it would bypass the
required freeze step. A corrupt record blocks reuse with a diagnostic; it is
not silently overwritten. Concurrent controllers cannot replace each other's
resolved choices. Preserve the settings directory with the run's provenance.

This is a development feature in the HybAs draft. Bactologix controls and
exact-commit execution compatibility still need release integration.
