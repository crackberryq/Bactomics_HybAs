"""End-to-end storage contracts with synthetic tool stubs only."""
import gzip
from pathlib import Path
from test_workflow import prepare, run, calls


def test_generated_streams_native_sidecars_and_prior_results_survive(tmp_path):
    _,env=prepare(tmp_path,{'sample':'hybrid'})
    run(tmp_path,env)
    isolate=tmp_path/'sample'
    raw=isolate/'work/inputs_merged/R1.merged.uncompressed.fastq'
    compressed=isolate/'work/inputs_merged/illumina_R1.fastq.gz'
    assert raw.is_file()
    # The uncompressed merge is retained independently of downstream transformations.
    assert raw.read_bytes().startswith(b'@fixture')
    filt=isolate/'work/inputs_filtered/ont_filtered.fastq.gz'
    assert Path(str(filt)+'.tmp.fastq').read_bytes()==gzip.decompress(filt.read_bytes())
    sam=isolate/'work/mapping/illumina.paired_qc.sam'
    assert sam.is_file() and sam.stat().st_size
    native=isolate/'work/assembly/unicycler'
    sidecar=native/'native-sentinel.bin';sidecar.write_bytes(b'old native intermediate')
    assembly=(native/'assembly.fasta').read_bytes()
    unicycler=[c for c in calls(tmp_path) if c['tool']=='unicycler'][-1]['args']
    assert unicycler[unicycler.index('--keep')+1]=='3'
    assert '--no-clear-after' in unicycler[unicycler.index('--spades_options')+1]
    run(tmp_path,env,'--forcerun','unicycler')
    matches=list((tmp_path/'retained').glob('*/files/sample/work/assembly/unicycler/native-sentinel.bin'))
    assert len(matches)==1 and matches[0].read_bytes()==b'old native intermediate'
    assert (matches[0].parent/'assembly.fasta').read_bytes()==assembly
    saved_sam=list((tmp_path/'retained').glob('*/files/sample/work/mapping/illumina.paired_qc.sam'))
    assert saved_sam and all(p.stat().st_size for p in saved_sam)
    before=set((tmp_path/'retained').iterdir())
    run(tmp_path,env)
    assert set((tmp_path/'retained').iterdir())==before
