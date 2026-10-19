# Contributing

Keep HybAs directly runnable through `Hybas.smk` and Snakemake-managed environments. Preserve the separation of AWARE target election, independent platform actions and confidence. A confidence change must not silently become a read-cleaning policy change.

Run `python -m pytest tests -q` before proposing a change. The frozen v8.8 controller fixture documents the default compatibility contract; do not edit it to make a changed policy pass. Add regression coverage for an actual defect or required behavior, and disclose biological validation limits.

Keep execution-affecting settings in rule parameters/dependencies, retain durable logs on failure, and avoid mutating another rule's inputs. New optional stages must record whether they executed and which artifact supplied the final output.

Do not commit workstation paths, sample data, database downloads, active run outputs or credentials. Documentation must describe the files actually shipped. The historical phase 1–3 scientific validation bundle is separate from this repository's software tests.
