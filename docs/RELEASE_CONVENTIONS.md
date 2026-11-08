# Bactomics update conventions

From 2026-10-29, new updates follow the naming and release-document structure
used by [BaSIS AWARE 0.4.0a2](https://github.com/crackberryq/Bactomics_BaSIS/releases/tag/v0.4.0a2).
This is the standing convention for future HybAs, BaSIS and Bactologix work.
It does not claim that every repository's release automation already enforces it.

## Names and versions

- Use the product's canonical name: **HybAs AWARE**, **BaSIS AWARE** or
  **Bactologix**. AWARE is an edition name, not a version suffix.
- Use `MAJOR.MINOR.PATCH` with `aN`, `bN` or `rcN` for alpha, beta or release
  candidates. Keep each product's existing version sequence independently.
- Name a release `<Product> <version> — <stage>`, for example
  `HybAs AWARE 10.0.0a1 — Alpha`. A stable release can omit the stage suffix.
- Use `v<version>` for the Git tag, for example `v10.0.0a1`.
- Mark development work explicitly as unreleased. A proposed version, branch
  or successful test run is not a published release.
- Use one consistent version across runtime metadata, README, citation,
  changelog, release notes and the release catalogue when preparing publication.
  A user-supplied legacy version label is not proof of the source version.

## Release documents and source layout

Keep one `docs/RELEASE_<version>.md` per new release, linked from the README.
Follow BaSIS's order: identity and status; concrete changes; update and test
instructions; verification and release scope. Record the exact tested source
commit and CI evidence in the published release. Include known gaps and any
compatibility or migration steps relevant to that update.

For standalone workflows, retain a top-level Snakefile, a pinned controller
environment, configuration/examples, workflow helpers, rule environments,
documentation, tests and CI. HybAs already separates these roles; its current
paths are documented in [REPOSITORY_STRUCTURE.md](REPOSITORY_STRUCTURE.md).
Move existing paths only with an explicit migration and reference updates.
Bactologix remains an installable application with its wheel/source artifacts;
the shared update format does not change its distribution type.

BaSIS's release line is `basis-aware`; the matching intended HybAs release-line
name is `hybas-aware`. Feature branches remain separate. The current HybAs work
stays in its existing draft branch until integration; this convention does not
create or redirect a release branch.

## Publication and compatibility

Publish only from the exact commit that passes the required release checks.
Keep alpha/beta/RC releases marked as prereleases. Published tags and historical
release identities remain fixed; another update gets another version.

Bactologix discovery and execution compatibility are separate: a published
release must also have a supported adapter/contract for its exact source.
Update that compatibility record as part of release integration. Preserve
existing installations, frozen run configurations, outputs and provenance.

The next HybAs target is documented in [RELEASE_10.0.0a1.md](RELEASE_10.0.0a1.md).
