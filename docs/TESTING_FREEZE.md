# HybAs 10 development feature freeze

The testing candidate includes declared Illumina/ONT/PacBio input routes, Auto settings, regular AWARE selection/cleaning, identification review, stop-after-AWARE and verified continuation, optional extended read preservation, portable derived samples, assignment-support statistics, branded reports and computational provenance. No additional analysis routes should enter this testing cycle; fixes should be recorded against the frozen pair.

## Acceptance matrix

| Area | Cases | Acceptance evidence |
| --- | --- | --- |
| Inputs | Illumina PE; repeated/complementary Illumina; supported mate-pair; ONT; Illumina+ONT; PacBio HiFi/CLR with and without Illumina | Correct declared route, chemistry, pairing, preprocessing and library identities; no silent omission |
| AWARE controls | Regular defaults; preservation on/off; statistics on/off; automatic/manual selection | Existing decisions unchanged when extensions are disabled; optional outputs match enabled settings |
| Stage boundaries | Identification pause; stop after AWARE; direct assembly; continue saved parent | No downstream completion claims before execution; verified upstream reuse |
| Derived samples | Primary/species selection; copied bundle; hybrid parent yielding single-platform child | Actual retained inputs determine route; parent identity retained; no accidental second filtering |
| Failure and recovery | Interrupted jobs; changed inputs/evidence; missing environments; insufficient storage | Clear failure reason, retained previous results, no misleading completed state |
| Reports and Methods | Desktop/mobile; identification-only; AWARE-only; assembled; failed continuation | Actual target distinguished from recommendation; stage-specific prose and downloadable evidence |
| Installation | Fresh WSL/Linux; existing workstation; isolated app workspace | Matching controller and tool environments; reproducible installation and resume |
| Reference-data qualification | Representative environmental bacterial datasets with known references; legacy v9 comparator | Predeclared comparisons of assembly correctness/completeness, read retention, resources and differences; inspect discordant results |

## Test record

Record the paired source commits, app version, machine/OS, environment identities, database versions/checksums, input identities, configuration, expected outcome, observed outcome, logs and report paths. Keep software fixture results separate from reference-data accuracy evaluations. A successful software test does not qualify biological accuracy.

Run each dataset in its own output directory. Compare regular AWARE and enabled extensions using the same frozen original inputs and settings. Do not treat overlapping primary/species groups or repeated computational runs as independent biological replicates. Define analytical acceptance thresholds before inspecting comparative outcomes.

## Packaging and release gates

The Bactologix repository provides `scripts/package_development_pair.py`. It builds source archives and an application wheel from explicit immutable commits, plus a pair manifest and checksums. It never publishes, changes the stable channel or grants managed execution approval.

Before managed GUI testing of a newly packaged HybAs source, complete its exact-commit native adapter CI and add the corresponding prerelease compatibility entry through the official release process. Before release qualification, complete the reference-data matrix and document unresolved limitations. Existing source-bound runs remain tied to their original versions.
