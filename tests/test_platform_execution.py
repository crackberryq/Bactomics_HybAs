"""Generic synthetic workflow contracts, not biological accuracy benchmarks."""
import gzip
import json
from pathlib import Path

import pytest
import yaml

from test_workflow import prepare, run, calls, settings_snapshots, computational_runs
from test_library_execution import manifest_run


def platform_run(tmp, platform, chemistry='hifi', short=False, polish='none', kraken=True):
    cfg, env = prepare(tmp, {'sample': 'illumina_only' if short else 'ont_only'}, polish=polish, kraken=kraken)
    libs = []
    if short:
        libs.append(dict(library_id='pe', platform='illumina', library_type='paired_end', orientation='fr',
            read_groups=[dict(read_group_id='pe1', run_id='run1',
                r1=str(tmp/'sample/raw/illumina/demo_R1.fastq.gz'), r2=str(tmp/'sample/raw/illumina/demo_R2.fastq.gz'))]))
    path = tmp/'long reads.fq.gz'
    path.write_bytes(gzip.compress(b'@duplicate\n'+b'A'*1500+b'\n+\n'+b'I'*1500+b'\n'))
    lib = dict(library_id='long', platform=platform,
               read_groups=[dict(read_group_id='long1', run_id='longrun', reads=str(path))])
    if platform == 'pacbio': lib['chemistry'] = chemistry
    libs.append(lib)
    mode = ('illumina_' + platform) if short else platform + '_only'
    manifest = tmp/'inputs.yaml'
    manifest.write_text(yaml.safe_dump(dict(schema_version=1, samples=[dict(sample_id='sample', mode=mode, libraries=libs)])))
    config = yaml.safe_load(cfg.read_text());config.update(input_manifest=str(manifest), sample_sheet='')
    cfg.write_text(yaml.safe_dump(config))
    return cfg, env, manifest


@pytest.mark.parametrize('chemistry,short,polish', [('hifi',False,'none'), ('clr',False,'tripolish'), ('hifi',True,'tripolish'), ('clr',True,'none')])
def test_pacbio_routes(tmp_path, chemistry, short, polish):
    _, env, _ = platform_run(tmp_path, 'pacbio', chemistry, short, polish)
    run(tmp_path, env)
    invoked = calls(tmp_path)
    flye = next(c for c in invoked if c['tool']=='flye')
    assert ('--pacbio-hifi' if chemistry=='hifi' else '--pacbio-raw') in flye['args']
    assert not any(c['tool'] in {'medaka_consensus','racon','filtlong','unicycler','spades.py'} for c in invoked)
    summary = json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())
    assert summary['assembler']=='flye' and summary['long_read_platform']=='pacbio'
    assert summary['pacbio']['chemistry']==chemistry
    assert summary['pacbio']['raw']['bases']==1500
    assert summary['pacbio_kraken_decision']=='passthrough'
    report = tmp_path/'sample/reports/hybas_research'
    assert 'PacBio filtering' in (report/'hybas_report.html').read_text()
    audit = json.loads((report/'aware_pacbio_filter_audit.json').read_text())
    assert audit['status'] == 'completed' and audit['integrity_ok']
    assert not (report/'aware_ont_filter_audit.json').exists()
    decision = (tmp_path/'sample/metadata/kraken_decision.tsv').read_text()
    assert 'pacbio_selected_target_percent\t100.0000' in decision
    assert '\nont_selected_target_percent\t' not in decision
    assert not (tmp_path/'sample/work/inputs_filtered/pacbio_filtered.fastq.gz').exists()
    assert summary['executed_stages']==(['polypolish'] if short and polish=='tripolish' else [])
    before = calls(tmp_path);run(tmp_path,env)
    assert calls(tmp_path)==before
    assert computational_runs(tmp_path)[-1]['status']=='no_work'


@pytest.mark.parametrize('short', [False, True])
def test_ont_manifest_routes(tmp_path, short):
    _, env, _ = platform_run(tmp_path, 'ont', short=short, polish='tripolish')
    run(tmp_path, env)
    summary=json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())
    assert summary['assembler']==('spades' if short else 'unicycler')
    assert summary['long_read_platform']=='ont'
    assert summary['executed_stages']==['racon1','racon2','medaka']+(['polypolish'] if short else [])


@pytest.mark.parametrize('kind,orientation,preprocessing', [
    ('paired_end','rf','pretrimmed'), ('paired_end','ff','pretrimmed'),
    ('mate_pair','rf','pretrimmed'), ('mate_pair','fr','pretrimmed'),
    ('mate_pair','rf','nextera_mate_pair')])
def test_declared_library_preparations(tmp_path,kind,orientation,preprocessing):
    _,env,manifest=manifest_run(tmp_path,polish='none')
    data=yaml.safe_load(manifest.read_text());lib=data['samples'][0]['libraries'][0]
    lib.update(library_type=kind,orientation=orientation,preprocessing=preprocessing)
    manifest.write_text(yaml.safe_dump(data));run(tmp_path,env)
    dataset=json.loads((tmp_path/'sample/trimmed/libraries/dataset.yaml').read_text())
    expected='paired-end' if kind=='paired_end' else 'mate-pairs'
    assert any(d['type']==expected and d['orientation']==orientation for d in dataset)
    assert sum(c['tool']=='fastp' for c in calls(tmp_path))==1
    assert any(c['tool']=='nxtrim' for c in calls(tmp_path))==(preprocessing=='nextera_mate_pair')


def test_pacbio_ont_polisher_override_rejected(tmp_path):
    _,env,_=platform_run(tmp_path,'pacbio',polish='medaka')
    result=run(tmp_path,env,ok=False)
    assert 'incompatible with PacBio' in result.stdout
    assert not calls(tmp_path) and not settings_snapshots(tmp_path)


def test_pacbio_failure_resume_preserves_preparation(tmp_path):
    _, env, _ = platform_run(tmp_path, 'pacbio', 'clr')
    failed = dict(env, HYBAS_TEST_FAIL_TOOL='flye')
    run(tmp_path, failed, ok=False)
    before = calls(tmp_path)
    assert computational_runs(tmp_path)[-1]['status'] == 'failed'
    run(tmp_path, env, '--rerun-incomplete')
    after = calls(tmp_path)
    assert sum(c['tool']=='kraken2' for c in before) == sum(c['tool']=='kraken2' for c in after)
    assert sum(c['tool']=='flye' for c in after) == 2
    assert (tmp_path/'sample/reports/hybas_research/hybas_report.html').is_file()
