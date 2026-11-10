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

This is not a guarantee that every transient created internally by every third-party executable survives. In particular, Dev's Flye 2.9.6 deletes some internal consensus/polishing alignments, temporary bubbles/consensus files and unfiltered contigs without a supported keep flag; its remaining native output directory is retained. A BUSCO keep-all internal-scratch guarantee has not been established. Streaming bytes that were never written to a file and kernel/runtime scratch are not captured. Reconstructible SQLite indexes and lock files used by HybAs bookkeeping can still be removed; the original reads, assignment ledgers and report evidence remain.

Hard termination, disk failure or explicit user deletion cannot be prevented by this adapter. Existing files already deleted by older releases cannot be recovered. Manual Snakemake cleanup/delete commands are explicit destructive operations, not ordinary execution, and are outside the rerun adapter.

## Development additions

The Dev branch preserves each replaced computational HTML report under `provenance/runtime/report_history/`, including completion-callback refreshes outside scheduled jobs. It additionally retains per-library cleaned inputs in `trimmed/libraries/cleaned_inputs/attempt-*`, failed AWARE staging as `reports/aware.attempt-*` with `failure.json`, and generated bootstrap selector reports under the AWARE bundle's `bootstrap_reports/`. A failed AWARE bundle is rejected for continuation. Failed Flye/SPAdes working directories have visible unique attempt names under `work/assembly/`.

Tests use synthetic products and tool-contract stubs to verify retention, failure/retry behavior and unchanged workflow routing. They do not establish scientific accuracy.

## Conda controller identity

The executing `snakemake.__version__` identifies the supported engine, matching the CLI. Conda's known `0.0.0` distribution-metadata placeholder (or absent distribution metadata) cannot replace that identity. An unsupported runtime or a concrete conflicting distribution version stops execution. New run version records also use the executing module version. Release publishing requires the native Conda controller retention tests and a synthetic end-to-end version-record check, in addition to the existing regression and browser checks.
