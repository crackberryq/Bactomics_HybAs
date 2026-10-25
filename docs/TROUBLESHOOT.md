# Troubleshooting

Read the failed rule and its `logs/<rule>.log` before changing inputs. Snakemake's exit-code summary alone may not contain the tool's underlying error.

| Symptom | Check |
|---|---|
| No isolates | Populate the header-only `samples.tsv`, or select the intended config |
| Unpaired Illumina lanes | Match every `_R1` filename with its `_R2` lane; remove duplicate compressed/uncompressed copies |
| Invalid taxid | Use a positive numeric ID or YAML null/empty; literal `None` is invalid |
| Missing database | Verify paths; run the explicit database admin target from the README |
| BUSCO directory exists but preflight fails | A directory/stamp is insufficient; dataset metadata, cutoff files and HMMs are required |
| SPAdes memory failure | Read durable assembly logs; check `assembly_mem_mb`, `spades_memory_gb`, host/WSL RAM and concurrent jobs |
| Severe swapping | Check Kraken database size, `kraken_mem_mb` and the profile total; keep `kraken_jobs: 1` unless capacity is measured |
| Lock after interruption | Ensure the old process has stopped, then `snakemake -s Hybas.smk --unlock` |
| Failed/incomplete output | Resume with the same command and `--rerun-incomplete`; keep valid upstream outputs |
| Broken Conda environment | Identify the exact environment from the error. Unlocking does not delete or repair a completed environment |
| Missing validation workflow | The historical phase 1–3 bundle is separate and not shipped here |

Do not clear `.snakemake`, Kraken outputs or all work directories as a general restart method. Report-only refreshes and polishing failures should preserve valid classifications.

For a damaged environment already marked complete, remove/rebuild only that identified environment after stopping active jobs. For a failed environment creation, Snakemake manages the incomplete deployment on rerun. Do not use `--cleanup-metadata` to pretend incomplete output is valid.

Synthetic tests in certain hosted sandboxes cannot read child-process statistics because `/proc` uses a different PID namespace. The README documents a test-copy-only benchmark switch; production benchmarking is unchanged.
