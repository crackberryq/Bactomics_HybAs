Computational performance and execution history
==============================================

Every execution records its own UTC timestamps, elapsed time, hardware and runner
details, requested resources, and job attempts. The tables distinguish measured
CPU time and memory from requested threads and memory. Archived JSON, TSV and raw
benchmarks are retained under ``provenance/runtime/computational_runs/``.

Elapsed run time includes planning and scheduling; parallel job durations must
not be interpreted as elapsed run time. The largest job RSS peak is not a
simultaneous workflow memory peak. Native benchmarks are sampled measurements;
short or failed jobs may have no CPU or memory sample. Missing values remain
explicit. A successful conditional rule can still have skipped its tool; consult
the run's existing skip and execution records for that distinction.
