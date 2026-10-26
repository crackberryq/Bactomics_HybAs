"""Real Snakemake tests with synthetic tool-contract stubs, plus failure/resume."""
import gzip
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml

ROOT=Path(__file__).resolve().parents[1]
SNAKE=shutil.which('snakemake') or str(Path(sys.executable).with_name('snakemake'))
TOOLS=['conda','pigz','fastqc','NanoPlot','fastp','filtlong','kraken2','extract_kraken_reads.py',
       'unicycler','minimap2','racon','medaka_consensus','bwa','samtools','polypolish','quast.py','seqkit','multiqc','prokka','busco']


def prepare(tmp, modes, polish='tripolish', pp=True, kraken=False):
    tmp.mkdir(parents=True,exist_ok=True)
    binpath=tmp/'bin';binpath.mkdir()
    for name in TOOLS:(binpath/name).symlink_to(ROOT/'tests/mock_tool.py')
    for name,mode in modes.items():
        files=[]
        if mode!='ont_only':files += ['illumina/demo_R1.fastq.gz','illumina/demo_R2.fastq.gz']
        if mode!='illumina_only':files += ['nanopore/demo.fastq.gz']
        for file in files:
            path=tmp/name/'raw'/file;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(gzip.compress(b'@fixture\n'+b'A'*1200+b'\n+\n'+b'I'*1200+b'\n',mtime=0))
    samples=tmp/'samples.tsv';samples.write_text('isolate\n'+'\n'.join(modes)+'\n')
    db=tmp/'db'/'kraken'
    if kraken:
        db.mkdir(parents=True)
        for name in ['hash.k2d','opts.k2d','taxo.k2d']:(db/name).write_text('fixture')
    config=dict(base_dir=str(tmp),sample_sheet=str(samples),threads=2,run_kraken=kraken,
        kraken_decision_mode='aware' if kraken else 'off',run_busco=False,run_prokka=False,run_multiqc=True,
        run_polypolish=pp,polish_mode=polish,racon_rounds=2,db_root='db',kraken_db_name='kraken',
        assembly_mem_mb=1024,spades_memory_gb=1,polishing_mem_mb=512)
    path=tmp/'config.yaml';path.write_text(yaml.safe_dump(config))
    env=dict(os.environ,PATH=str(binpath)+os.pathsep+str(Path(SNAKE).parent)+os.pathsep+os.environ['PATH'],
             HYBAS_TEST_CALLS=str(tmp/'calls.jsonl'))
    return path,env


def run(tmp, env, *extra, ok=True):
    snakefile=ROOT/'Hybas.smk'
    if env.get('HYBAS_TEST_DISABLE_BENCHMARKS') == '1':
        # Some hosted execution sandboxes expose a /proc from another PID
        # namespace. Only the test copy omits resource sampling there; CI uses
        # the shipped Snakefile, including benchmarks, without this setting.
        copy=tmp/'workflow';copy.mkdir(exist_ok=True)
        for item in ROOT.iterdir():
            if item.name != 'Hybas.smk' and not (copy/item.name).exists():
                (copy/item.name).symlink_to(item)
        snakefile=copy/'Hybas.smk'
        content=re.sub(r'^    benchmark:[^\n]*\n(?:        [^\n]+\n)?', '', (ROOT/'Hybas.smk').read_text(), flags=re.M)
        if not snakefile.exists() or snakefile.read_text()!=content:snakefile.write_text(content)
    command=[SNAKE,'-s',str(snakefile),'--configfile',str(tmp/'config.yaml'),'--cores','2','--latency-wait','1',*extra]
    result=subprocess.run(command,cwd=tmp,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
    (tmp/'last-run.log').write_text(result.stdout)
    if ok:assert result.returncode==0,result.stdout[-18000:]
    else:assert result.returncode!=0,result.stdout
    return result


def calls(tmp):
    path=tmp/'calls.jsonl'
    return [json.loads(x) for x in path.read_text().splitlines()] if path.exists() else []


def computational_runs(tmp):
    return [json.loads(path.read_text()) for path in sorted(
        (tmp/'provenance/runtime/computational_runs').glob('*/run.json'))]


@pytest.mark.parametrize('mode,polish,pp,expected,stages',[
    ('hybrid','tripolish',True,'CCGT',['racon1','racon2','medaka','polypolish']),
    ('hybrid','tripolish',False,'CCG',['racon1','racon2','medaka']),
    ('illumina_only','tripolish',True,'AT',['polypolish']),
    ('illumina_only','none',True,'AA',[]),
    ('ont_only','tripolish',True,'CCG',['racon1','racon2','medaka']),
    ('hybrid','racon',True,'ACC',['racon1','racon2']),
])
def test_end_to_end_routes(tmp_path,mode,polish,pp,expected,stages):
    _,env=prepare(tmp_path,{'sample':mode},polish,pp)
    run(tmp_path,env)
    text=(tmp_path/'sample/work/final/assembly.final.fasta').read_text().strip()
    assert text.endswith(expected)
    summary=json.loads((tmp_path/'sample/metadata/run_summary.json').read_text())
    computational=computational_runs(tmp_path)
    assert len(computational)==1
    assert computational[0]['status']=='success'
    assert computational[0]['collection_errors']==[]
    jobs=computational[0]['jobs']
    assert all(job['status']=='success' for job in jobs)
    assert any(job['rule']=='unicycler' for job in jobs)
    if env.get('HYBAS_TEST_DISABLE_BENCHMARKS')!='1':
        measured=[job for job in jobs if job['benchmark_path']]
        assert measured
        assert all(job['measurements'] for job in measured)
        archive=tmp_path/'provenance/runtime/computational_runs'/computational[0]['run_id']/measured[0]['benchmark_archive']
        assert 'input_size_mb' in archive.read_text().splitlines()[0]  # Extended even with custom profiles.
    assert summary['executed_stages']==stages
    assert float(summary['quast_n50'])>0
    assert summary['final_source_stage']==(stages[-1] if stages else 'unicycler')
    if 'polypolish' not in stages:
        assert not any(c['tool']=='bwa' for c in calls(tmp_path))
        assert not (tmp_path/'provenance/conda/polishing.conda_list.txt').exists()
    if 'medaka' not in stages:assert not (tmp_path/'provenance/conda/medaka.conda_list.txt').exists()
    # A no-op rerun must be a genuine no-op, including admin/provenance rules.
    before=len(calls(tmp_path));result=run(tmp_path,env)
    assert len(calls(tmp_path))==before,result.stdout
    reruns=computational_runs(tmp_path)
    assert len(reruns)==2
    assert reruns[0]==computational[0]
    assert reruns[1]['status']=='no_work'
    assert reruns[1]['jobs']==[]


def test_batch_failure_resume_and_report_refresh(tmp_path):
    _,env=prepare(tmp_path,{'good':'hybrid','bad':'hybrid'},kraken=True)
    failing=dict(env,HYBAS_TEST_FAIL_TOOL='polypolish',HYBAS_TEST_FAIL_ISOLATE='/bad/')
    run(tmp_path,failing,'--keep-going',ok=False)
    failed=computational_runs(tmp_path)[0]
    assert failed['status']=='failed'
    assert any(job['rule']=='polypolish' and job['status']=='failed' for job in failed['jobs'])
    assert all(not job['measurements'] for job in failed['jobs'] if job['status']=='failed')
    good=tmp_path/'good/work/final/assembly.final.fasta'
    assert good.exists()
    original=good.stat().st_mtime_ns
    classifications=[c for c in calls(tmp_path) if c['tool']=='kraken2']
    run(tmp_path,env)
    resumed=computational_runs(tmp_path)
    assert resumed[0]==failed
    assert resumed[1]['status']=='success'
    assert not any(job['rule'].startswith('kraken2_') for job in resumed[1]['jobs'])
    assert good.stat().st_mtime_ns==original
    assert [c for c in calls(tmp_path) if c['tool']=='kraken2']==classifications
    before=calls(tmp_path)
    run(tmp_path,env,'--forcerun','run_summary')
    new=calls(tmp_path)[len(before):]
    assert {c['tool'] for c in new} <= {'multiqc'},new


def test_configuration_change_reuses_classification(tmp_path):
    cfg,env=prepare(tmp_path,{'sample':'hybrid'},kraken=True)
    run(tmp_path,env)
    before=calls(tmp_path)
    data=yaml.safe_load(cfg.read_text());data['target_taxid']=12;cfg.write_text(yaml.safe_dump(data))
    run(tmp_path,env)
    new=calls(tmp_path)[len(before):]
    assert not any(c['tool']=='kraken2' for c in new)
    decision=(tmp_path/'sample/metadata/kraken_decision.tsv').read_text()
    assert 'selected_target_taxid\t12\n' in decision
    assert any(c['tool']=='extract_kraken_reads.py' for c in new)


def test_empty_mapping_is_explicit_skip(tmp_path):
    _,env=prepare(tmp_path,{'sample':'illumina_only'})
    env['HYBAS_TEST_NO_ALIGNMENTS']='1'
    run(tmp_path,env)
    report=json.loads((tmp_path/'sample/metadata/polish_execution.json').read_text())
    assert report['executed_stages']==[]
    assert report['skipped_stages']['polypolish']=='no_mapped_alignments'
    assert report['final_source_stage']=='unicycler'


@pytest.mark.parametrize('mode',['hybrid','ont_only','illumina_only'])
@pytest.mark.parametrize('polish',['none','racon','medaka','tripolish'])
def test_dag_modes(tmp_path,mode,polish):
    _,env=prepare(tmp_path,{'sample':mode},polish)
    result=run(tmp_path,env,'-n')
    assert 'MissingInputException' not in result.stdout
    assert computational_runs(tmp_path)==[]


def test_computational_quiet_progress_and_report(tmp_path):
    _,env=prepare(tmp_path,{'sample':'illumina_only'},polish='none',pp=False)
    run(tmp_path,env,'--quiet','progress')
    before=computational_runs(tmp_path)
    assert before[0]['status']=='success'
    assert all(job['status']=='success' for job in before[0]['jobs'])
    assert 'Computational record' in (tmp_path/'provenance/runtime/computational_report.html').read_text()
    run(tmp_path,env,'--report',str(tmp_path/'report.zip'))
    assert computational_runs(tmp_path)==before  # Report rendering is not an execution run.
    from zipfile import ZipFile
    with ZipFile(tmp_path/'report.zip') as archive:
        html=[archive.read(name).decode() for name in archive.namelist() if name.endswith('.html')]
    assert any('Computational record' in item for item in html)
    run(tmp_path,env,'--touch','--forcerun','run_summary')
    assert computational_runs(tmp_path)==before


def test_invalid_pair_blocks_entire_batch(tmp_path):
    _,env=prepare(tmp_path,{'good':'hybrid','bad':'hybrid'})
    (tmp_path/'bad/raw/illumina/demo_R2.fastq.gz').unlink()
    result=run(tmp_path,env,ok=False)
    assert 'Unpaired Illumina lanes' in result.stdout
    assert calls(tmp_path)==[]


def test_setup_without_samples(tmp_path):
    _,env=prepare(tmp_path,{})
    run(tmp_path,env,'setup_envs','--config','setup_only=true')
    assert not any(c['tool'] in {'unicycler','kraken2','polypolish'} for c in calls(tmp_path))


def test_oversized_reservation_blocks_before_processing(tmp_path):
    cfg,env=prepare(tmp_path,{'sample':'hybrid'},kraken=True)
    data=yaml.safe_load(cfg.read_text());data['kraken_mem_mb']=20000;cfg.write_text(yaml.safe_dump(data))
    result=run(tmp_path,env,'-n',ok=False)
    assert 'scheduling budget' in result.stdout
    assert calls(tmp_path)==[]
