"""Scheduler/control-flow acceptance with stubs, never scientific accuracy tests."""
import json
import os
from pathlib import Path
import yaml
from test_workflow import prepare, run, calls


def test_busco_content_change_invalidates_even_with_identical_stat(tmp_path):
    config,env=prepare(tmp_path,{'sample':'illumina_only'},polish='none',pp=False)
    dataset=tmp_path/'db/busco/bacteria_odb10';(dataset/'hmms').mkdir(parents=True)
    for name in ['dataset.cfg','scores_cutoff','lengths_cutoff','hmms/fixture.hmm']:(dataset/name).write_text('AAAA')
    cfg=yaml.safe_load(config.read_text());cfg['run_busco']=True;config.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    first=json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())['busco_dataset_sha256']
    before=calls(tmp_path);run(tmp_path,env);assert calls(tmp_path)==before
    hmm=dataset/'hmms/fixture.hmm';st=hmm.stat();hmm.write_text('BBBB');os.utime(hmm,ns=(st.st_atime_ns,st.st_mtime_ns))
    run(tmp_path,env)
    fresh=calls(tmp_path)[len(before):]
    assert sum(x['tool']=='busco' for x in fresh)==1
    assert not any(x['tool']=='unicycler' for x in fresh)
    after=json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())['busco_dataset_sha256']
    assert after!=first


def test_medaka_upgrade_preserves_legacy_consensus_and_new_receipt_matches(tmp_path):
    _,env=prepare(tmp_path,{'sample':'ont_only'},polish='medaka',pp=False)
    old=tmp_path/'sample/work/polishing/medaka/consensus.fasta';old.parent.mkdir(parents=True)
    old.write_text('>legacy\nAAAA\n')
    run(tmp_path,env)
    assert old.read_text()=='>legacy\nAAAA\n'
    receipt=json.loads((old.parent/'execution.json').read_text())
    summary=json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())
    assert summary['medaka_model_effective']=='fixture-model.tar.gz'
    assert summary['medaka_model_sha256']==receipt['model']['sha256']
    assert '/current/consensus.fasta' in json.loads((tmp_path/'sample/metadata/polish_execution.json').read_text())['final_source_path']
    assert Path(receipt['attempt'],'native/consensus.fasta').is_file()
    receipt['consensus_sha256']='0'*64
    (old.parent/'execution.json').write_text(json.dumps(receipt))
    rejected=run(tmp_path,env,str(tmp_path/'sample/metadata/polish_execution.json'),
                 '--forcerun','polish_execution',ok=False)
    assert 'Medaka execution receipt does not match' in rejected.stdout
    assert old.read_text()=='>legacy\nAAAA\n'
