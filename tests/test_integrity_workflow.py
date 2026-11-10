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


def test_v9_original_model_tags_drive_auto_selection_and_mixed_rerun_stops(tmp_path):
    import gzip
    config,env=prepare(tmp_path,{'sample':'ont_only'},polish='medaka',pp=False)
    source=tmp_path/'sample/raw/nanopore/demo.fastq.gz'
    record=b'@fixture RG:Z:run_modelA_barcode01 model_version_id=modelA\n'+b'A'*1200+b'\n+\n'+b'I'*1200+b'\n'
    source.write_bytes(gzip.compress(record))
    cfg=yaml.safe_load(config.read_text());cfg['medaka_model']='';config.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    merged=tmp_path/'sample/work/inputs_merged/ont_merged.fastq.gz'
    with gzip.open(merged,'rb') as f:assert f.read()==record
    receipt=json.loads((tmp_path/'sample/work/polishing/medaka/execution.json').read_text())
    assert receipt['model']['selection']=='automatic_from_complete_read_metadata'
    assert receipt['model']['input_metadata']['basecaller_models']==['modelA']
    previous_calls=calls(tmp_path)
    source.write_bytes(gzip.compress(record+record.replace(b'@fixture ',b'@other ').replace(b'modelA',b'modelB')))
    run(tmp_path,env,ok=False)
    medaka_log=(tmp_path/'sample/logs/medaka.log').read_text()
    assert 'Automatic Medaka selection requires one model identity' in medaka_log
    assert not any(x['tool']=='medaka_consensus' for x in calls(tmp_path)[len(previous_calls):])
    assert Path(receipt['attempt'],'native/consensus.fasta').is_file()
