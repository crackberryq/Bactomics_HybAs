"""Stage boundaries and unchanged upstream reuse, using software-only fixtures."""
import json
from pathlib import Path
import sys
import yaml
import pytest

from test_workflow import prepare, run, calls, ROOT
sys.path.insert(0, str(ROOT/'workflow/lib'))
import identification_review as review


def stage(tmp_path, mode='hybrid'):
    path, env = prepare(tmp_path, {'sample': mode}, polish='none', pp=False, kraken=True)
    cfg=yaml.safe_load(path.read_text())
    cfg.update(kraken_review_before_filtering=True, kraken_decision_mode='manual',
               require_target_taxid_for_kraken=True, target_taxid='')
    path.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    p=tmp_path/'sample/reports/identification/review.json'
    return env,p,json.loads(p.read_text())


def decision(tmp_path,p,doc,mode='manual'):
    target=doc['recommendation']['taxid']
    value={'schema':review.DECISION_SCHEMA,'mode':mode,'taxid':target,
           'evidence_path':str(p),'evidence_sha256':review.digest(p)}
    path=tmp_path/'decision.json';path.write_text(json.dumps(value))
    return ['--config','kraken_review_decision='+str(path),'kraken_review_decision_sha256='+review.digest(path)]


@pytest.mark.parametrize('mode', ['hybrid','illumina_only','ont_only'])
def test_identify_review_continue_reuses_kraken(tmp_path,mode):
    env,p,doc=stage(tmp_path,mode)
    assert doc['schema']==review.SCHEMA
    assert doc['recommendation']['taxid']
    assert not (tmp_path/'sample/work/final/assembly.final.fasta').exists()
    initial=calls(tmp_path)
    assert not any(x['tool'] in {'unicycler','extract_kraken_reads.py','fastp','filtlong'} for x in initial)
    first_count=sum(x['tool']=='kraken2' for x in initial)
    extra=decision(tmp_path,p,doc)
    run(tmp_path,env,*extra)
    assert (tmp_path/'sample/work/final/assembly.final.fasta').exists()
    assert sum(x['tool']=='kraken2' for x in calls(tmp_path))==first_count
    assert json.loads(p.read_text())==doc
    before=len(calls(tmp_path));run(tmp_path,env,*extra)
    assert len(calls(tmp_path))==before


def test_review_cannot_be_bypassed_by_explicit_assembly_target(tmp_path):
    env,p,doc=stage(tmp_path)
    result=run(tmp_path,env,str(tmp_path/'sample/work/final/assembly.final.fasta'),ok=False)
    assert 'Awaiting researcher review' in result.stdout


def test_changed_evidence_and_database_block_continuation(tmp_path):
    env,p,doc=stage(tmp_path)
    extra=decision(tmp_path,p,doc)
    raw=Path(doc['files'][0]['path']);raw.write_text(raw.read_text()+'\n')
    result=run(tmp_path,env,*extra,ok=False)
    assert 'Saved Kraken evidence changed' in result.stdout


def test_identification_does_not_require_assembly_database_or_manual_taxid(tmp_path):
    path,env=prepare(tmp_path,{'sample':'illumina_only'},polish='none',kraken=True)
    cfg=yaml.safe_load(path.read_text());cfg.update(kraken_review_before_filtering=True,run_busco=True,
        kraken_decision_mode='manual',require_target_taxid_for_kraken=True,target_taxid='')
    path.write_text(yaml.safe_dump(cfg));run(tmp_path,env)
    assert (tmp_path/'sample/reports/identification/report.html').is_file()


def test_auto_continuation_clears_original_manual_override(tmp_path):
    env,p,doc=stage(tmp_path,'illumina_only')
    run(tmp_path,env,*decision(tmp_path,p,doc,'aware'))
    fields=dict(line.split('\t',1) for line in (tmp_path/'sample/metadata/kraken_decision.tsv').read_text().splitlines()[1:])
    assert fields['kraken_decision_mode']=='aware'
    assert fields['manual_target_taxid']==''
    assert fields['selected_target_taxid']==doc['recommendation']['taxid']


def test_database_change_invalidates_review_context(tmp_path):
    env,p,doc=stage(tmp_path)
    extra=decision(tmp_path,p,doc)
    (tmp_path/'db/kraken/taxo.k2d').write_text('a changed software fixture')
    result=run(tmp_path,env,*extra,ok=False)
    assert 'settings or source changed' in result.stdout


@pytest.mark.parametrize('platform,chemistry,short', [
    ('pacbio','hifi',False), ('pacbio','clr',False),
    ('pacbio','hifi',True), ('pacbio','clr',True), ('ont','hifi',True)])
def test_manifest_review_platform_and_frozen_resume(tmp_path, platform, chemistry, short):
    from test_platform_execution import platform_run
    from test_workflow import settings_snapshots
    cfg,env,_=platform_run(tmp_path,platform,chemistry,short)
    data=yaml.safe_load(cfg.read_text());data.update(kraken_review_before_filtering=True)
    cfg.write_text(yaml.safe_dump(data))
    run(tmp_path,env)
    p=tmp_path/'sample/reports/identification/review.json';doc=json.loads(p.read_text())
    assert doc['platforms']=={'illumina':short,platform:True}
    assert platform+'_pct' in doc['candidates'][0]
    assert ('PacBio %' if platform=='pacbio' else 'ONT %') in p.with_name('report.html').read_text()
    initial=calls(tmp_path)
    assert not any(x['tool'] in {'flye','spades.py','unicycler','extract_kraken_reads.py','fastp','filtlong'} for x in initial)
    snapshots={str(p):p.read_bytes() for p in settings_snapshots(tmp_path)}
    blocked=run(tmp_path,env,str(tmp_path/'sample/work/final/assembly.final.fasta'),ok=False)
    assert 'Awaiting researcher review' in blocked.stdout
    run(tmp_path,env,*decision(tmp_path,p,doc,'aware'))
    assert (tmp_path/'sample/work/final/assembly.final.fasta').is_file()
    assert sum(x['tool']=='kraken2' for x in calls(tmp_path))==sum(x['tool']=='kraken2' for x in initial)
    assert {str(p):p.read_bytes() for p in settings_snapshots(tmp_path)}==snapshots


def test_dual_library_review_preserves_library_identity(tmp_path):
    from test_library_execution import manifest_run
    cfg,env,_=manifest_run(tmp_path,polish='none',kraken=True)
    data=yaml.safe_load(cfg.read_text());data['kraken_review_before_filtering']=True
    cfg.write_text(yaml.safe_dump(data));run(tmp_path,env)
    p=tmp_path/'sample/reports/identification/review.json';doc=json.loads(p.read_text())
    initial=calls(tmp_path)
    assert not any(x['tool'] in {'spades.py','fastp'} for x in initial)
    run(tmp_path,env,*decision(tmp_path,p,doc,'aware'))
    assert sum(x['tool']=='kraken2' for x in calls(tmp_path))==sum(x['tool']=='kraken2' for x in initial)
    dataset=json.loads((tmp_path/'sample/trimmed/libraries/dataset.yaml').read_text())
    assert len(dataset)==2
