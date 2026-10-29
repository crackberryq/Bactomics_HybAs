# Computational data

Collection is part of normal HybAs v10 execution. There is no analysis config
toggle to disable it. Native benchmarking covers every concrete rule; aggregate
targets such as `all` perform no computation and have no benchmark. Extended
benchmarks are enabled in the workflow as well as the default profile.

## Files and history

Each controller execution writes a unique directory under
`<base_dir>/provenance/runtime/computational_runs/<UTC-time>-<id>/`:

- `run.json`: run status, UTC start/end, monotonic elapsed seconds, controller
  hardware/runner information, source/config hashes, requested resources and all
  job attempts with measured resource values.
- `events.jsonl`: append-only job start/end records, flushed as events arrive.
- `resolved_settings.json`: the settings request, selected values/reasons and
  sample plans for this invocation. `run.json` links the immutable settings
  identity and includes its decisions alongside current scheduler limits.
- `jobs.tsv`: one row per benchmark repetition or unmeasured attempt, exported
  when execution ends. Missing measurements are `NA` (JSON uses `null`).
- `benchmarks/`: copies of fresh native benchmark files, archived when a job
  succeeds, before a later invocation can overwrite its usual benchmark path.

`provenance/runtime/computational_report.html` presents the run history and stage
measurements. It is also registered in the native Snakemake report. Existing
per-sample `benchmarks/*.tsv` remain available; global/provenance rules use
`provenance/runtime/benchmarks/`.

Each run has a collapsible resolved-settings section. The shared checksummed
snapshot under `provenance/runtime/resolved_settings/` is reused for an unchanged
request; the per-invocation copy preserves its association with this run. See
[Auto settings and resume](AUTO_SETTINGS.md) for precedence and re-detection.

Successful, failed and resumed invocations get separate records. Retries retain
separate attempts. A no-work invocation records zero new job attempts and never
imports the earlier run's costs. Quiet progress output is supported. Report
rendering, dry runs and touch operations do not create execution measurements.

## Measurements and interpretation

| Field | Meaning |
| --- | --- |
| `elapsed_seconds` | Controller execution elapsed time, including planning, scheduling and waiting; not a sum of concurrent job durations. |
| `wall_seconds` | Native benchmark duration for that job/repetition. |
| `cpu_seconds` | Native sampled user + system CPU time. |
| `max_rss_mib`, `max_vms_mib`, `max_uss_mib`, `max_pss_mib` | Native sampled memory peaks, in MiB. |
| `io_read_mib`, `io_write_mib` | Native sampled I/O counters. The pinned 9.27.1 implementation converts bytes to MiB. These are not output-directory sizes. |
| `mean_cpu_load` | Snakemake's native `mean_load` value, retained without reinterpretation. |
| `input_size_mib` | Input size from the extended benchmark, when available. |
| `threads_requested`, `resources_requested` | Scheduler allocations; these are not measured CPU/RAM consumption. |
| `controller_elapsed_seconds` | Time between controller job announcements, including execution overhead; null if completion is only reconciled at workflow end. |
| `measured_job_cpu_seconds_sum` | Sum of available measurements only. Coverage counts indicate missing attempts. |
| `largest_job_peak_rss_mib` | Largest individual job RSS peak, not simultaneous total workflow peak memory. |

Hardware information describes the **controller** host: OS/architecture, Python,
Snakemake, logical/physical CPU counts, accessible CPU affinity, OS-reported RAM
and initial free space. Scheduler resource limits are recorded separately.
Remote worker hardware, GPU utilization and energy use are not measured here.

Snakemake resource measurements are sampled, not exact accounting. Short jobs,
unavailable OS counters and terminated processes can leave partial measurements.
Failed jobs generally have no native benchmark: their controller timing/status is
retained, while unavailable CPU and memory values stay null. Pure Python `run:`
rules have the native benchmark's process-level measurement limitations. Job
success means the rule finished; a conditional rule may have skipped its tool.
Use existing execution/skip metadata for that distinction.

## Failure boundaries

Collection begins after Snakemake announces execution and selects its executor,
after the Snakefile has loaded. Syntax/configuration errors that prevent loading have no execution record.
Normal success and failure hooks finalize records; exceptional shutdown leaves an
explicit `incomplete` status. A hard kill or power loss can leave `running` records
with no end time, which must not be treated as successful. Job starts and previous
completions are retained as written. No software can finalize after a hard kill.

Collection I/O errors produce a warning and are included in the record when the
record can still be saved. They do not silently become zero resource usage.
History is not automatically deleted or rewritten on resume.

Tests use synthetic tool stubs. Hosted containers with incompatible PID namespaces
can strip native benchmarks from temporary test copies via the documented
`HYBAS_TEST_DISABLE_BENCHMARKS` test setting; production code has no such bypass.
CI exercises the normal workflow and checks archived extended benchmark headers.
