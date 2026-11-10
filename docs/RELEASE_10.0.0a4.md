# HybAs 10.0.0a4

This development patch allows an explicit additional swap allowance to reach the Snakemake scheduler. Previously, the native settings resolver capped the entire scheduler budget to physical RAM even when Bactologix supplied a larger RAM-plus-swap budget.

The frozen `bactologix_compute.swap_mb` setting declares the additional allowance. Snakemake's `--resources mem_mb` must include both a positive RAM portion and this allowance. The resolver subtracts swap before resolving native physical-RAM settings, preserves its existing 80% Auto policy and physical-RAM ceilings, then adds the allowance only to the effective scheduler budget. Assembly, SPAdes, polishing and other tool settings are not inflated by swap. Kraken's database-size reservation is not reduced to fit.

The native effective budget may be lower than the requested Bactologix budget because the native physical-RAM policy remains more conservative. The controller log shows physical RAM, additional swap, and the effective scheduler budget separately. Resolved-settings provenance records the split; historical records are preserved on resume. A changed allowance receives a separate settings identity.

The allowance is a scheduling declaration, not an operating-system limit or a measurement of memory consumed. Bactologix checks currently usable RAM, free Linux swap and container limits before launch. Standalone users must verify that capacity themselves: this native resolver does not verify or allocate swap. Paging can substantially slow execution, and a successful dry run does not establish runtime performance or analytical accuracy.

Use a new run with this source. Software regression coverage includes insufficient and sufficient scheduler budgets, malformed policies, unchanged tool RAM settings, resume provenance, and real Snakemake 9.27.1 dry runs with inert sparse fixtures. No analysis tools run in the new scheduling tests, and no assembly, filtering, polishing or quality thresholds change.
