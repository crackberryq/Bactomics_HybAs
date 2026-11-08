"""Synthetic software fixtures for AWARE endpoints; no biological validation."""
import json
import sys
from pathlib import Path
import pytest
import yaml
from test_workflow import prepare, run, calls, ROOT
sys.path.insert(0, str(ROOT/'workflow/lib'))
from identification_review import digest
from aware_bundle import verify_bundle, CONTINUATION_SCHEMA
from aware_statistics import wilson, resample

@pytest.mark.parametrize('mode', ['illumina_only','ont_only','hybrid'])
def test_stop_and_explicit_continue(tmp_path, mode):
    path,env=prepare(tmp_path, {'sample':mode}, polish='none', pp=False, kraken=True)
    cfg=yaml.safe_load(path.read_text());cfg['aware_run_until']='aware';path.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    manifest=tmp_path/'sample/reports/aware/manifest.json'
    doc=verify_bundle(manifest)
    assert doc['groups'][0]['assembly_eligible']
    assert not any(c['tool'] in {'unicycler','fastp','filtlong'} for c in calls(tmp_path))
    assert not (tmp_path/'sample/work/final/assembly.final.fasta').exists()
    result=run(tmp_path,env,str(tmp_path/'sample/work/final/assembly.final.fasta'),ok=False)
    assert 'Stop after AWARE' in result.stdout
    before=sum(c['tool']=='kraken2' for c in calls(tmp_path))
    continuation=tmp_path/'continue.json';continuation.write_text(json.dumps(dict(schema=CONTINUATION_SCHEMA,bundle_path=str(manifest),bundle_sha256=digest(manifest))))
    run(tmp_path,env,'--config','aware_assembly_continuation='+str(continuation),'aware_assembly_continuation_sha256='+digest(continuation))
    assert (tmp_path/'sample/work/final/assembly.final.fasta').exists()
    assert sum(c['tool']=='kraken2' for c in calls(tmp_path))==before
    assert digest(manifest)==json.loads(continuation.read_text())['bundle_sha256']


def test_wilson_and_bootstrap_counts():
    import random
    assert wilson(0,0) is None
    assert wilson(50,100)==pytest.approx([.4038315,.5961685])
    assert wilson(0,100)[1]==pytest.approx(.0369935)
    rng=random.Random(1729)
    draws=[resample(rng,{'a':40,'b':60}) for _ in range(3000)]
    assert all(sum(d.values())==100 for d in draws)
    assert sum(d['a'] for d in draws)/len(draws)==pytest.approx(40,abs=.3)


@pytest.mark.parametrize('platform,chemistry,short', [('pacbio','hifi',True),('pacbio','clr',False),('ont','hifi',True)])
def test_manifest_extended_endpoint(tmp_path, platform, chemistry, short):
    from test_platform_execution import platform_run
    path,env,_=platform_run(tmp_path, platform, chemistry, short)
    cfg=yaml.safe_load(path.read_text());cfg.update(aware_run_until='aware',aware_preserve_groups=True,aware_statistics=True)
    path.write_text(yaml.safe_dump(cfg));run(tmp_path,env)
    doc=verify_bundle(tmp_path/'sample/reports/aware/manifest.json')
    primary=doc['groups'][0]
    assert primary['assembly_eligible']
    libs=primary['manifest']['samples'][0]['libraries']
    long=next(l for l in libs if l['platform']==platform)
    if platform=='pacbio':assert long['chemistry']==chemistry
    assert doc['statistics']['stability']['same_target_fraction']==1
    assert not any(c['tool'] in {'flye','spades.py','unicycler'} for c in calls(tmp_path))


def test_review_then_stop(tmp_path):
    from test_identification_review import decision
    path,env=prepare(tmp_path,{'sample':'illumina_only'},polish='none',pp=False,kraken=True)
    cfg=yaml.safe_load(path.read_text());cfg.update(aware_run_until='aware',kraken_review_before_filtering=True)
    path.write_text(yaml.safe_dump(cfg));run(tmp_path,env)
    review=tmp_path/'sample/reports/identification/review.json'
    assert not (tmp_path/'sample/reports/aware/manifest.json').exists()
    run(tmp_path,env,*decision(tmp_path,review,json.loads(review.read_text()),'aware'))
    assert (tmp_path/'sample/reports/aware/manifest.json').exists()
    assert not (tmp_path/'sample/work/final/assembly.final.fasta').exists()


def test_retained_complementary_library_metadata(tmp_path):
    from aware_bundle import LibraryWriter
    from aware_preservation import Writers
    from input_manifest import validate_manifest
    plan={'mode':'auto','selection':{'selected_mode':'dual_illumina_complementary'},'libraries':[
        dict(library_id='pe',platform='illumina',library_type='paired_end',orientation='fr',insert_size_bp=400,read_groups=[dict(run_id='run_pe',read_group_id='lane_pe')]),
        dict(library_id='mp',platform='illumina',library_type='mate_pair',orientation='rf',insert_size_bp=3000,preprocessing='nextera_mate_pair',read_groups=[dict(run_id='run_mp',read_group_id='lane_mp')])]}
    mapper=LibraryWriter(tmp_path,plan);writers=Writers(tmp_path)
    try:
        for i in range(2):
            name=f'hybasL{i}G0R1'
            record=(name,('@'+name+'\nACGT\n+\nIIII\n').encode(),4)
            mapper.write(writers,'primary','illumina',name,[record,record])
    finally:writers.close()
    result=mapper.manifest('primary','child')
    assert validate_manifest(result,base_dir=tmp_path)['selections'][0]['selected_mode']=='dual_illumina_complementary'
    assert result['samples'][0]['libraries'][1]['preprocessing']=='nextera_mate_pair'
    assert result['samples'][0]['libraries'][1]['read_groups'][0]['run_id']=='run_mp'

    long_id=mapper.manifest('primary','s'*100+'_aware_primary')['samples'][0]['sample_id']
    assert len(long_id)<=100 and long_id.startswith('s')


def test_completed_bundle_records_actual_filtering(tmp_path):
    path,env=prepare(tmp_path,{'sample':'illumina_only'},polish='none',pp=False,kraken=True)
    cfg=yaml.safe_load(path.read_text());cfg.update(aware_run_until='aware',aware_statistics=True)
    path.write_text(yaml.safe_dump(cfg));run(tmp_path,env)
    doc=verify_bundle(tmp_path/'sample/reports/aware/manifest.json')
    assert doc['filtering']['selected_target_taxid']
    html=(tmp_path/'sample/reports/aware/report.html').read_text()
    assert 'Selected target and actual filtering' in html and '95% Wilson interval' in html
