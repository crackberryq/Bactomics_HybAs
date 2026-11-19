# HybAs AWARE 10.0.0a7 — Alpha

This release fixes the BUSCO database preparation step. BUSCO 6 can extract a
dataset such as `bacteria_odb10` under `lineages/bacteria_odb10` in its download
directory. HybAs previously checked the directory one level above the extracted
dataset and reported it as incomplete.

HybAs now locates the requested lineage and verifies its configuration, HMM
profiles and cutoff files before marking the database ready. In Bactologix Dev,
select this HybAs version and run **Download and prepare** on the BUSCO lineage
card. Choose **Use downloaded** after the preparation job succeeds.
