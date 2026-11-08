"""Software routing contracts, not analytical validation."""
import json
import yaml
import pytest
from test_workflow import prepare, run, calls, ROOT
from test_controller import controller
from aware_modes import resolve

@pytest.mark.parametrize('mode', ['off', 'taxonomy', 'auto', 'split'])
def test_mode_execution(tmp_path, mode):
    path,env=prepare(tmp_path, {'sample':'illumina_only'}, polish='none', pp=False, kraken=True)
    cfg=yaml.safe_load(path.read_text());cfg['aware_mode']=mode;path.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    invoked=[c['tool'] for c in calls(tmp_path)]
    assert ('kraken2' in invoked)==(mode!='off')
    assert 'extract_kraken_reads.py' not in invoked
    assert (tmp_path/'sample/work/final/assembly.final.fasta').exists()==(mode!='split')
    if mode!='off':
        doc=json.loads((tmp_path/'sample/reports/aware/manifest.json').read_text())
        assert doc['filtering']['aware_mode']==mode
        assert doc['filtering']['decision_confidence']=='high'
        assert doc['filtering']['platform_agreement_reason']=='single_platform_not_assessed'
        assert 'Experimental' in (tmp_path/'sample/reports/aware/report.html').read_text()


def test_ambiguous_auto_review_and_resume(tmp_path):
    from test_identification_review import decision
    path,env=prepare(tmp_path, {'sample':'illumina_only'}, polish='none', pp=False, kraken=True)
    stub=tmp_path/'bin/kraken2';stub.unlink()
    stub.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
p=sys.argv[sys.argv.index('--report')+1]
Path(p).write_text('100\\t2\\t0\\tR\\t1\\troot\\n50\\t1\\t1\\tG\\t11\\tExample\\n50\\t1\\t1\\tG\\t22\\tOther\\n')
print('C\\tfixture\\t11\\t1200|1200\\t11:4 |:| 11:4')
''');stub.chmod(0o755)
    cfg=yaml.safe_load(path.read_text());cfg['aware_mode']='auto';path.write_text(yaml.safe_dump(cfg))
    run(tmp_path,env)
    report=tmp_path/'sample/reports/identification/review.json';doc=json.loads(report.read_text())
    assert doc['awaiting_review']
    assert not any(c['tool'] in {'unicycler','fastp','extract_kraken_reads.py'} for c in calls(tmp_path))
    result=run(tmp_path,env,str(tmp_path/'sample/work/final/assembly.final.fasta'),ok=False)
    assert 'Awaiting researcher review' in result.stdout
    run(tmp_path,env,*decision(tmp_path,report,doc,'aware'))
    assert (tmp_path/'sample/work/final/assembly.final.fasta').exists()


def test_reconstruction_unavailable_and_legacy_unchanged():
    with pytest.raises(ValueError,match='not implemented'):resolve({'aware_mode':'reconstruct'})
    cfg={'kraken_decision_mode':'manual','target_taxid':'11'}
    assert resolve(cfg)==cfg


def test_single_platform_is_not_disagreement(tmp_path):
    c=controller(ROOT/'Hybas.smk');c['AWARE_MODE']='auto'
    left,right=tmp_path/'l',tmp_path/'r'
    left.write_text('99\t99\t99\tS\t11\tExample\n');right.write_text('')
    d=c['compute_taxonomic_confidence'](left,right,'11')
    assert d['decision_confidence']=='high'
    assert d['platform_agreement_reason']=='single_platform_not_assessed'


def test_ancestor_is_not_a_competing_population(tmp_path):
    c=controller(ROOT/'Hybas.smk');c['AWARE_MODE']='taxonomy'
    left,right=tmp_path/'l',tmp_path/'r'
    left.write_text('100\t100\t1\tG\t10\tExample\n99\t99\t99\tS\t11\tExample one\n1\t1\t1\tS\t12\tExample two\n');right.write_text('')
    d=c['compute_taxonomic_confidence'](left,right,'11')
    assert d['decision_confidence']=='high'
    assert d['runner_up_taxid']=='12'
