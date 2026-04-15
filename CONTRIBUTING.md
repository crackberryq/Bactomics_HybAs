## Contributing to HybAs Aware

Thank you for your interest in improving **HybAs Aware**.

At this stage, the repository is maintained as a **controlled development project with limited outside contributions**. External feedback is welcome, especially when it improves reproducibility, documentation quality, workflow clarity, and bug reporting discipline.

This document explains how contributions should be proposed, discussed, and prepared before code is merged.

---

## 1. Contribution model

HybAs Aware is not currently managed as an open-ended community project with unrestricted drive-by changes.

The preferred contribution model is:

* open an issue first
* discuss major changes before opening a pull request
* keep changes scoped and technically justified
* prioritize reproducibility and repository consistency over large speculative refactors

Bug reports, reproducibility issues, documentation corrections, and well-scoped improvements are welcome.

---

## 2. Before opening an issue or pull request

Please review the existing repository documentation first.

Relevant files include:

* `README.md`
* `INSTALL.md`
* `USAGE.md`
* `CONFIGURATION.md`
* `TROUBLESHOOTING.md`
* `REPOSITORY_STRUCTURE.md`

If your question is already answered there, update your local setup first and re-check the issue before submitting a report.

---

## 3. What kinds of contributions are most useful

The most helpful contributions at this stage are:

### Bug reports

Especially when they include:

* exact command used
* execution context
* relevant config values
* isolate layout summary
* first failing rule
* full or relevant error output
* expected behavior versus observed behavior

### Reproducibility reports

Examples:

* environment resolution differences
* platform-specific failures
* version mismatch problems
* `.condarc` or channel-related behavior differences
* unexpected divergence between CLI and Streamlit behavior

### Documentation improvements

Examples:

* unclear path explanations
* inconsistent terminology
* missing caveats
* confusing command examples
* places where the docs imply behavior more strongly than the code supports

### Small, targeted workflow improvements

These are most useful when they:

* preserve existing terminology
* preserve the isolate-centric layout model
* do not silently change public workflow behavior
* include documentation updates where needed

---

## 4. What to discuss before making major changes

Please open an issue before proposing changes that affect any of the following:

* repository structure
* isolate layout assumptions
* controller environment model
* `.condarc` execution behavior
* run mode logic
* polishing logic
* taxon-aware decision logic
* downstream phase behavior
* output naming
* provenance or metadata schema
* public-facing workflow terminology

These are core parts of the repository contract and should not be changed casually.

---

## 5. Expectations for bug reports

A good bug report should help reproduce the issue without guesswork.

Please include as many of the following as possible:

### Environment context

* operating system
* shell environment
* whether you are using the pinned Ubuntu baseline or another platform
* Conda version
* Snakemake version
* Python version
* whether the controller environment was created from `environment.yaml`

### Execution context

* whether you ran via CLI or Streamlit
* working directory used
* whether `CONDARC` was exported from the repository path
* exact command used

### Workflow context

* whether the issue occurred in the main workflow or the downstream workflow
* isolate name
* sample-sheet row used
* whether the run was standard mode, admin mode, or CI mode

### Error context

* first failing rule
* first real error message
* relevant log excerpt
* whether the run was interrupted previously

Please focus on the **first real failure**, not only on later downstream symptoms.

---

## 6. Expectations for pull requests

If you open a pull request, please keep it structured and limited in scope.

### Pull requests should:

* address a clearly defined problem
* avoid mixing unrelated changes
* preserve current path conventions unless explicitly discussed first
* keep wording consistent with repository terminology
* update documentation when behavior changes
* avoid silent changes to expected outputs or config semantics

### Pull requests should not:

* rename major workflow files without prior discussion
* reorganize the isolate layout unilaterally
* move the `validation/` workflow code just because of folder naming preference
* repurpose `environment.yaml` into a full rule-tool environment file
* remove `.condarc`-based execution assumptions without an agreed migration plan

---

## 7. Repository terminology that should remain consistent

Please preserve the current public terminology unless a change is explicitly approved.

Use:

* **HybAs Aware** for the public project name
* **HybAs** as the short form where context is clear
* **Downstream Computational Genome Analysis Workflow** as the public name for the workflow implemented under `validation/`
* **controller environment** for the environment created from `environment.yaml`
* **rule environments** for the per-rule Conda environments created by Snakemake

Avoid introducing alternate names for the same components inside documentation or pull requests unless the rename is intentional and discussed first.

---

## 8. Repository structure assumptions contributors should respect

The current repository model assumes:

* `bactomics/` is the top-level workspace
* isolates are direct children of `bactomics/`
* `bactomics/hybas/` contains workflow code and control assets
* `bactomics/hybas/validation/` contains downstream workflow code
* `bactomics/<isolate>/validation/` contains downstream generated outputs

Contributions should respect the distinction between:

* workflow code paths
* isolate input/output paths

Do not collapse these concepts together in code or documentation.

---

## 9. Configuration and behavior changes

Configuration changes should be made carefully.

If you change:

* `config.yaml`
* `validation/validation_config.yaml`
* `ci/config.ci.yaml`
* sample sheet assumptions
* environment pins

then you should also review whether corresponding updates are needed in:

* `README.md`
* `INSTALL.md`
* `USAGE.md`
* `CONFIGURATION.md`
* `TROUBLESHOOTING.md`

Configuration drift between code and docs is one of the easiest ways to degrade repository quality.

---

## 10. Documentation changes

Documentation contributions are welcome and strongly encouraged when they improve clarity without inventing unsupported behavior.

Good documentation changes:

* clarify current behavior based on code evidence
* remove ambiguity
* tighten path descriptions
* correct command examples
* flag caveats where users are likely to make mistakes

Documentation changes should not:

* invent features not supported by the repository
* claim compatibility not yet verified
* add benchmark claims without evidence
* add installation steps that are not actually part of the repository model

---

## 11. Reproducibility and evidence expectations

When reporting or proposing technical changes, prefer evidence over impression.

Helpful evidence includes:

* exact file paths
* config snippets
* sample-sheet rows
* environment details
* first failing command or rule
* before/after behavior

When possible, use repository terminology consistently and point to the relevant file rather than describing behavior loosely.

---

## 12. Style expectations for contributions

### Keep changes focused

A single pull request should solve one coherent problem where possible.

### Preserve naming stability

Do not rename files, folders, config keys, or outputs without prior discussion.

### Avoid speculative cleanup

Large stylistic refactors are less helpful than narrowly justified improvements.

### Keep docs and code aligned

If a change affects behavior, update the documentation in the same pull request.

### Prefer explicitness

Path changes, config changes, and workflow-behavior changes should be obvious from the diff.

---

## 13. Support boundaries

HybAs Aware is currently documented as a controlled-development repository.

That means:

* outside contributions may be reviewed selectively
* not every feature request will be accepted
* maintainers may prioritize reproducibility and scope discipline over expansion

This is normal for a repository that is still being stabilized for public technical use.

---

## 14. Security and sensitive data

When opening issues or pull requests:

* remove private sample identifiers unless they are safe to share
* sanitize local absolute paths where possible
* do not commit sensitive data
* do not commit large raw sequencing data
* avoid exposing internal-only infrastructure details unless necessary

For documentation examples, prefer sanitized relative paths.

---

## 15. Minimal contribution workflow

A good default process is:

1. open an issue
2. describe the problem clearly
3. wait for discussion if the change is structural or behavioral
4. prepare a focused branch
5. update documentation if needed
6. open a pull request with a clear summary

This is especially important for changes that affect workflow logic, repository structure, or public-facing terminology.

---

## 16. Good issue topics at this stage

Examples of useful issues include:

* incorrect path assumptions in docs
* broken command examples
* mismatch between config behavior and README wording
* reproducibility differences between CLI and Streamlit runs
* `.condarc`-related execution problems
* CI profile confusion
* missing caveats for downstream workflow prerequisites
* version string inconsistencies across files

---

## 17. Contributions that need extra caution

Please discuss first before proposing:

* major dependency changes
* environment model changes
* moving or renaming `Hybas.smk`
* changing isolate path assumptions
* changing output directory conventions
* changing metadata schema
* changing workflow version semantics
* replacing or heavily restructuring the downstream workflow phases

These changes can ripple across the entire repository and documentation package.

---

## 18. Status of this contribution policy

This `CONTRIBUTING.md` reflects the current project stance:

* controlled development
* limited outside contributions
* bug reports and reproducibility issues welcome
* documentation corrections welcome
* major changes should be discussed before a pull request

It is intentionally lightweight, but it is still strict about preserving repository consistency and evidence-based documentation.
