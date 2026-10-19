# Configuration

`config.yaml` contains portable defaults. Local paths are relative to the launch working directory; `base_dir` identifies the parent of all isolate folders. `samples.tsv` supplies isolate IDs and optional taxid, BUSCO lineage and Medaka model overrides.

| Setting | Meaning |
|---|---|
| `threads` | Requested rule threads, bounded by Snakemake `--cores` |
| `keep_percent` | Filtlong retention setting |
| `racon_rounds` | Number of Racon rounds; zero omits Racon |
| `polish_mode` | `none`, `racon`, `medaka` or `tripolish`; see the README route table |
| `run_polypolish` | Enable the compatible Tripolish short-read stage |
| `run_busco`, `run_prokka`, `run_multiqc` | Optional module toggles |
| `kraken_decision_mode` | `off`, `manual`, `auto` or `aware`; authoritative classification policy |
| `run_kraken` | Legacy flag; a contradictory value produces a warning |
| `target_taxid` | Positive numeric manual override or empty/YAML null |
| `kraken_auto_passthrough_min_percent` | Per-platform pass-through threshold, default 95 |
| `kraken_auto_prefer_species_min_percent` | Species preference threshold, default 90 |
| `kraken_auto_use_both_platforms` | Cross-platform aggregation option, default true |
| `kraken_include_children` | Include descendant taxids during extraction |
| `kraken_aware_high_target_min_percent` | High confidence target support, default 95 |
| `kraken_aware_moderate_target_min_percent` | Moderate confidence support, default 80 |
| `kraken_aware_delta_min_percent` | High confidence margin, default 20 |
| `kraken_aware_moderate_delta_min_percent` | Moderate confidence margin, default 10 |
| `assembly_mem_mb` | Assembly memory reservation, default 8192 MB |
| `spades_memory_gb` | Forwarded SPAdes memory limit, default 8 GB |
| `polishing_mem_mb` | Polishing/mapping reservation, default 4096 MB |
| `kraken_mem_mb` | Optional explicit classification reservation; otherwise estimated from database size |
| `min_fasta_length` | Unicycler export threshold, default 1000 |
| `db_root`, `kraken_db_name` | Database locations beneath `base_dir` |
| `busco_download_lineages` | Optional list for explicit database preparation |
| `setup_only`, `admin_mode` | Preflight bypasses for the named setup/admin targets |
| `workflow_version` | Configuration label; source version is separately authoritative |

`profiles/default/config.yaml` controls total scheduling resources and restart/log defaults. Its example memory budget is 16000 MB, with one Kraken job at a time. Adjust it for the machine; resource declarations do not enforce operating-system limits.

Default AWARE election, platform-specific actions and conservative confidence remain separate. Confidence parameters affect confidence labels, not the cleaning action thresholds. `auto` is a Kraken mode; this release does not add a new polishing mode.

The shipped `examples/ci/config.ci.yaml` is a small input profile, not a biological validation benchmark. Software tests create isolated synthetic data and configuration dynamically.
