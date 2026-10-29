# Configuration

`config.example.yaml` contains portable defaults; copy it to `config.yaml`.
Supported run controls default to Auto; see [AUTO_SETTINGS.md](AUTO_SETTINGS.md)
for precedence, hardware detection and frozen resume behavior. Local paths are relative to the launch working directory; `base_dir` identifies the parent of all isolate folders. `samples.tsv` supplies isolate IDs and optional taxid, BUSCO lineage and Medaka model overrides.

| Setting | Meaning |
|---|---|
| `run_mode` | `auto` or an explicit compatible `hybrid`, `illumina_only`, `ont_only` mode |
| `threads` | `auto` or requested rule threads, bounded by available CPUs and Snakemake `--cores` |
| `total_mem_mb` | `auto` (80% of detected capacity) or an explicit scheduling budget |
| `settings_revision` | Change explicitly to resolve another plan without deleting earlier records |
| `keep_percent` | Filtlong retention setting |
| `racon_rounds` | Number of Racon rounds; zero omits Racon |
| `polish_mode` | `auto` (existing Tripolish preset), `none`, `racon`, `medaka` or `tripolish` |
| `run_polypolish` | `auto`/on/off for the compatible Tripolish short-read stage |
| `run_busco`, `run_prokka`, `run_multiqc` | `auto`/on/off; Auto retains enabled defaults and required-database checks |
| `kraken_decision_mode` | `off`, `manual`, `auto` or `aware`; authoritative classification policy |
| `run_kraken` | `auto` respects explicit `off` policy; conflicting explicit flags retain validation |
| `target_taxid` | Positive numeric manual override or empty/YAML null |
| `kraken_auto_passthrough_min_percent` | Per-platform pass-through threshold, default 95 |
| `kraken_auto_prefer_species_min_percent` | Species preference threshold, default 90 |
| `kraken_auto_use_both_platforms` | Cross-platform aggregation option, default true |
| `kraken_include_children` | Include descendant taxids during extraction |
| `kraken_aware_high_target_min_percent` | High confidence target support, default 95 |
| `kraken_aware_moderate_target_min_percent` | Moderate confidence support, default 80 |
| `kraken_aware_delta_min_percent` | High confidence margin, default 20 |
| `kraken_aware_moderate_delta_min_percent` | Moderate confidence margin, default 10 |
| `assembly_mem_mb` | `auto`: up to 8192 MB within the budget; or explicit MB |
| `spades_memory_gb` | `auto`: fit assembly reservation; or explicit GB |
| `polishing_mem_mb` | `auto`: up to 4096 MB within the budget; or explicit MB |
| `kraken_mem_mb` | `auto`: estimate from database size; or explicit MB |
| `min_fasta_length` | Unicycler export threshold, default 1000 |
| `db_root`, `kraken_db_name` | Database locations beneath `base_dir` |
| `busco_download_lineages` | Optional list for explicit database preparation |
| `setup_only`, `admin_mode` | Preflight bypasses for the named setup/admin targets |
| `workflow_version` | Configuration label; source version is separately authoritative |

`profiles/default/config.yaml` requests all available cores, with one Kraken job at a time. HybAs bounds scheduling using detected CPU/memory limits and any explicit Snakemake ceilings. Resource declarations do not enforce operating-system limits.

Default AWARE election, platform-specific actions and conservative confidence remain separate. Confidence parameters affect confidence labels, not the cleaning action thresholds. `kraken_decision_mode: auto` keeps its original algorithm meaning. `polish_mode: auto` selects the existing preset and input-compatible stages; it adds no new polishing algorithm.

The shipped `examples/ci/config.ci.yaml` is a small input profile, not a biological validation benchmark. Software tests create isolated synthetic data and configuration dynamically.
