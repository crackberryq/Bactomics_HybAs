# Generated-output retention

HybAs retains its generated products for later inspection. Current assemblies and per-stage polishing FASTAs stay in `<base_dir>/<isolate>/work/`. This update also retains generated uncompressed merged/extracted FASTQs, Filtlong's uncompressed output and the paired mapping SAM (`work/mapping/illumina.paired_qc.sam`). These are additional files; the existing compressed products and final paths remain unchanged.

Before a scheduled job replaces existing products, HybAs moves the old products into:

`<base_dir>/retained/<UTC timestamp>-<rule>-<unique ID>/files/<original relative path>`

The adjacent `manifest.json` records the original paths, archive paths, rule, wildcards and preservation state. It does not infer that an earlier attempt completed successfully. Keep this manifest with the files when copying an isolate's history. Global retained history is outside individual isolate folders, so copying only the isolate folder omits that history; copy the entire run/base directory for a complete run transfer.

This happens **before Snakemake deletes outputs**, including `directory()` outputs, on ordinary standalone runs and on retries within the same invocation. Declared outputs, logs and benchmarks are covered. Explicitly registered native tool directories (including Unicycler, QUAST, BUSCO, Prokka and MultiQC) are preserved with their sidecars. A dry run or an up-to-date job creates no archive. Failed outputs remain available at their active paths until a retry archives them. Failed or partial products must not be treated as completed results.

The adapter requires the supported pinned Snakemake versions (9.6.2 or 9.27.1), enables `keep-incomplete` and `notemp` independently of the chosen profile, and refuses to proceed if preservation fails. Archives use same-filesystem renames, never hardlinks: replacing a current file cannot change its archived bytes. An interrupted move leaves its manifest and any already-moved files for recovery. Output directories with symlinked ancestors, products outside the registered run/database roots, remote products and in-place input/output updates are refused when replacement would be required. Individual symlinks are retained as links with their original targets recorded; linked target contents are not duplicated. Relative links may no longer resolve after a move: this retains link metadata, not an independently portable copy of the target contents. Keep active outputs and their retention root on the same filesystem. Hardlinked regular files, including files inside output directories, are refused before any movement because surviving aliases could modify archived bytes.

The configured `DB_ROOT` is registered separately for database-administration outputs; its overwritten bookkeeping products are archived under `DB_ROOT/retained/`. Unregistered external paths remain refused. Partial `pip freeze` output is retained in uniquely named `.partial.*` files alongside the environment manifests.

## Native tool behavior

- Unicycler 0.5.1 runs with `--keep 3`, retaining its supported native intermediates, including SPAdes output and alignment files.
- SPAdes receives `--no-clear-after`, including when invoked through Unicycler. This prevents its supported post-stage cleanup.
- Prokka 1.15.6 and QUAST 5.3.0 use their documented `--debug` retention modes.
- Medaka keeps unique attempt directories and failed publication scratch. A failed Polypolish output is retained instead of removed by a shell trap.

These switches affect retention/diagnostics; scientific thresholds, tool selection, polishing rounds and input data are unchanged. No new per-stage BUSCO, annotation comparison or polishing-effect analysis is introduced by this storage update.

Retention increases disk use, particularly for uncompressed reads, alignments and native assembly intermediates. HybAs does not automatically purge these products. Move or delete selected old runs only through a deliberate storage-management action after checking your transfer.

## Limits

Capture covers the reviewed Python-managed cleanup and overwrite paths in the pinned native packages described below. It does not intercept every filesystem operation inside compiled executables. For example, Flye's bundled SAMtools removes internal sort spill chunks; those files are outside this adapter. Unwritten streams, kernel/runtime scratch and in-memory states are not captured. Reconstructible SQLite indexes and lock files used by HybAs bookkeeping can still be removed; the original reads, assignment ledgers and report evidence remain.

Hard termination, disk failure or explicit user deletion cannot be prevented by this adapter. Existing files already deleted by older releases cannot be recovered. Manual Snakemake cleanup/delete commands are explicit destructive operations, not ordinary execution, and are outside the rerun adapter.

## Development additions

The Dev branch preserves each replaced computational HTML report under `provenance/runtime/report_history/`, including completion-callback refreshes outside scheduled jobs. It additionally retains per-library cleaned inputs in `trimmed/libraries/cleaned_inputs/attempt-*`, failed AWARE staging as `reports/aware.attempt-*` with `failure.json`, and generated bootstrap selector reports under the AWARE bundle's `bootstrap_reports/`. A failed AWARE bundle is rejected for continuation. Failed Flye/SPAdes working directories have visible unique attempt names under `work/assembly/`.

Tests use synthetic products and tool-contract stubs to verify retention, failure/retry behavior and unchanged workflow routing. They do not establish scientific accuracy.

## Conda controller identity

The executing `snakemake.__version__` identifies the supported engine, matching the CLI. Conda's known `0.0.0` distribution-metadata placeholder (or absent distribution metadata) cannot replace that identity. An unsupported runtime or a concrete conflicting distribution version stops execution. New run version records also use the executing module version. Release publishing requires the native Conda controller retention tests and a synthetic end-to-end version-record check, in addition to the existing regression and browser checks.


## Native cleanup archives

New runs use a source-checked retention launcher for BUSCO 6.0.0; Dev also uses one for Flye 2.9.6. The launchers preserve generated files immediately before the reviewed native Python code removes or truncates them, then permit the original cleanup operation. This preserves the tools' normal control flow and arguments. Installed tool packages are not modified. An unsupported package version or changed reviewed source stops the launch instead of silently proceeding without retention.

Each native output root contains a visible `native-retained/` directory with separate sessions and generations. Keep this directory with the entire native output bundle. Session receipts identify the tool and reviewed source; event manifests retain the original relative path and archived location. These paths are included when the surrounding output directory is retained on a later workflow retry. Failed and partial artifacts remain evidence of an attempt, not a completed result.

BUSCO archives its reviewed Python-managed scratch removals and replaced directories. It retains the previous baseline of each shared native-job log before its first overwrite, then saves a separate stdout/stderr record for every native job, including captured output from failed jobs. These separate records preserve diagnostics that parallel jobs can overwrite in the original shared log files. The source-checked adapter is propagated to BUSCO's spawned workers through a private package copy stored with the session.

Dev's Flye archives the consensus alignment, polishing-iteration alignments, bubbles and consensus files, graph-edge alignment and unfiltered polished contigs removed by its reviewed Python cleanup. Earlier `minimap.stderr` and `params.json` revisions are retained before replacement. Existing alignment sidecars remain in the original native bundle. Retention does not add polishing rounds or per-stage quality calculations.

The native adapters refuse destructive operations outside their owned output root, symlinked output paths and hardlinked files. Archives contain independent previous bytes; no hardlinks are used. A persistent failure record makes retention failure fatal even if upstream code catches a cleanup exception. Disk-full or interrupted archives must be inspected rather than treated as a successful completed retention session.

Required release checks run against the pinned native package source, exercising cleanup and overwrite behavior with synthetic products and stubbed expensive computations. These checks test storage and control-flow contracts, not scientific accuracy. Existing runs frozen to an older source keep that source; use a new run with the updated pipeline for this retention behavior. Files already deleted by older runs cannot be recovered.
