# HybAs 10.0.0a11 development release

Open the HybAs research report for a single view of the chosen polishing strategy. Stage measurements and evidence files are also available through the trajectory and validation reports.

- Follow complete BUSCO markers, predicted protein lengths and assembly measurements in stage order. Projection views expose individual dimensions of the trajectory.
- Compare each retained stage directly with the final assembly. The edit-density figure places projected edits by transition and annotated rRNA windows on final contig coordinates. Hotspot details show read depth at assessed positions.
- Inspect short-protein counts, sequence complexity, variation in window counts and descriptive motif enrichment. Tables identify window boundaries and units. Adjacent-stage edit counts describe sequence changes. Their biological accuracy requires independent evidence.
- Reopen an older run in a compatible Bactologix Dev build to regenerate a report edition. When measurements are missing, the preview identifies the inputs and calculations required.
- Workflow setup includes the polishing-validation environment for eligible runs. Its deployment check exercises the alignment tools, Python libraries and plotting backend.

Retain the original output directory when comparing runs or following it into BaSIS. A fresh HybAs run has its own directory and source record.

Reports display the version of the renderer separately from the version that produced the run.
