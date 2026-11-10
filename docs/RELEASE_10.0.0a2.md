# HybAs AWARE 10.0.0a2

Alpha development update on `development`, prerelease tag `v10.0.0a2`. This patch preserves generated job products before Snakemake replacement, including same-invocation retries, and retains supported native intermediates, failed outputs, per-library cleaned reads and AWARE evidence. See [retention locations and limits](OUTPUT_RETENTION.md).

Existing assembly/filtering/polishing choices remain unchanged. The previous development capabilities and their scientific qualification limits still apply; this update does not provide the separate historical phase 1–3 analysis bundle. Native Flye internal cleanup remains an explicit retention limitation.

Install/update through the checked Bactologix source or this prerelease. Preserve existing runs and allow additional disk capacity. Validation uses real Snakemake storage/retry tests and synthetic workflow regressions; these do not establish biological accuracy or beta readiness.
