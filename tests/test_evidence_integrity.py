"""Synthetic integrity regressions; these do not measure biological accuracy."""
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/lib'))
from medaka_stage import execute, model_metadata
from resource_identity import inventory, verify
from workflow_support import quast_metric
from test_controller import controller


def fq(path, models):
    path.write_text(''.join(f'@read{i}' + (f' model_version_id={model}' if model else '') + '\nACGT\n+\nIIII\n' for i,model in enumerate(models)))


class NativeStub:
    def __init__(self, model):
        self.model = model
        self.calls = []
        self.fail = False
    def __call__(self, args, **kwargs):
        self.calls.append(args)
        if args[0] == 'medaka':
            return SimpleNamespace(stdout=str(self.model)+'\n')
        if self.fail:
            raise subprocess.CalledProcessError(17, args)
        native = Path(args[args.index('-o')+1])
        assert not native.exists(), 'must not receive old native cache'
        native.mkdir()
        (native/'consensus.fasta').write_text('>fixture\nACGT\n')
        return SimpleNamespace(returncode=0)


def test_medaka_scheduled_runs_never_reuse_cache_and_record_effective_model(tmp_path):
    reads,draft,model = [tmp_path/n for n in ('reads.fq','draft.fa','model.tar.gz')]
    fq(reads, ['modelA']); draft.write_text('>draft\nACGT\n'); model.write_text('model-v1')
    out, receipt = tmp_path/'medaka/consensus.fasta',tmp_path/'medaka/execution.json'
    out.parent.mkdir(); (out.parent/'calls_to_draft.bam').write_text('stale bam')
    (out.parent/'consensus_probs.hdf').write_text('stale hdf')
    stub = NativeStub(model)
    first = execute(reads,draft,out,receipt,2,runner=stub,known_models=())
    fq(reads,['modelB']); draft.write_text('>draft\nAAAA\n'); model.write_text('model-v2')
    second = execute(reads,draft,out,receipt,2,requested='reviewed-model',runner=stub,known_models=())
    assert first['attempt'] != second['attempt']
    assert Path(first['attempt'],'native/consensus.fasta').exists()
    assert first['inputs_sha256'] != second['inputs_sha256']
    assert first['model']['sha256'] != second['model']['sha256']
    assert second['model']['selection'] == 'explicit'
    assert second['model']['effective_path'] == str(model)
    assert json.loads(receipt.read_text()) == second
    assert all('-m' in args for args in stub.calls if args[0]=='medaka_consensus')
    prior = out.read_bytes(); stub.fail=True
    with pytest.raises(subprocess.CalledProcessError):
        execute(reads,draft,out,receipt,2,requested='reviewed-model',runner=stub,known_models=())
    assert out.read_bytes()==prior and json.loads(receipt.read_text())==second
    assert any(json.loads(p.read_text())['status']=='failed' for p in out.parent.glob('attempts/*/execution.json'))


@pytest.mark.parametrize('models', [[None], ['modelA',None], ['modelA']*100+['modelB']])
def test_auto_selection_rejects_missing_or_mixed_metadata_without_fallback(tmp_path,models):
    reads,draft = tmp_path/'reads.fq',tmp_path/'draft.fa'; fq(reads,models);draft.write_text('>d\nACGT\n')
    stub=NativeStub(tmp_path/'model')
    with pytest.raises(ValueError,match='requires one model identity'):
        execute(reads,draft,tmp_path/'out/consensus.fa',tmp_path/'out/execution.json',2,runner=stub,known_models=())
    assert not stub.calls



def test_busco_identity_detects_same_size_timestamp_preserving_replacement(tmp_path):
    dataset=tmp_path/'lineage';dataset.mkdir();hmm=dataset/'hmm';hmm.write_text('AAAA')
    (dataset/'.complete').write_text('unchanged marker')
    before=inventory(dataset);st=hmm.stat();hmm.write_text('BBBB');os.utime(hmm,ns=(st.st_atime_ns,st.st_mtime_ns))
    after=inventory(dataset)
    assert before['sha256']!=after['sha256']
    with pytest.raises(ValueError,match='changed'):
        verify(dataset,before['sha256'])
    receipt=tmp_path/'receipt.json';verify(dataset,after['sha256'],receipt)
    assert json.loads(receipt.read_text())==after
    (dataset/'.complete').write_text('new marker')
    assert inventory(dataset)['sha256']==after['sha256']


def test_legacy_confidence_ignores_ancestor_but_keeps_real_same_rank_rival(tmp_path):
    c=controller(ROOT/'Hybas.smk');c['AWARE_MODE']='legacy'
    left,right=tmp_path/'left',tmp_path/'right'
    for path in (left,right):path.write_text('100\t100\t0\tG\t10\tExample\n99\t99\t99\tS\t11\tExample one\n1\t1\t1\tS\t12\tExample two\n')
    result=c['compute_taxonomic_confidence'](left,right,'11')
    assert result['decision_confidence']=='high' and result['runner_up_taxid']=='12'
    for path in (left,right):path.write_text('100\t100\t0\tG\t10\tExample\n50\t50\t50\tS\t11\tExample one\n50\t50\t50\tS\t12\tExample two\n')
    assert c['compute_taxonomic_confidence'](left,right,'11')['awareness_state']=='taxonomically_ambiguous'


@pytest.mark.parametrize('complete,expected', [('0','PASS'),('98','PASS'),('101','WARN'),('NaN','WARN'),('not-known','WARN')])
def test_qc_recording_is_never_analytical_acceptance(tmp_path,complete,expected):
    content=(ROOT/'Hybas.smk').read_text().split('rule qc_assessment:',1)[1].split('\nrule multiqc_custom_content:',1)[0]
    block=textwrap.dedent(content.split('    run:\n',1)[1])
    values={}
    for name in ['run_summary','quast','busco','cov','kr_illumina','kr_ont','fastp_html','seq_illumina','seq_ont','seq_asm']:
        values[name]=str(tmp_path/name);Path(values[name]).write_text('fixture\n')
    Path(values['run_summary']).write_text('metric\tvalue\nassembly_bp\t100\nassembly_contig_count_fasta\t1\n')
    Path(values['quast']).write_text('Assembly\tfixture\nN50\t100\n')
    Path(values['busco']).write_text(f'C:{complete}%[S:0%,D:0%],F:0%,M:100%,n:100\n')
    Path(values['cov']).write_text('metric\tvalue\nont_coverage_x\t20\n')
    output=tmp_path/'qc.tsv';assembly=tmp_path/'assembly.fa';assembly.write_text('>x\nAAAA\n')
    scope=dict(os=os,json=json,input=SimpleNamespace(**values),output=[str(output)],wildcards=SimpleNamespace(isolate='sample'),
               metadata_dir=lambda _:str(tmp_path),final_asm=lambda _:str(assembly),PACBIO_SAMPLES=[],RUN_BUSCO=True,
               KRAKEN_PROFILE_ENABLED=True,RUN_MODE={'sample':'hybrid'},illumina_r1_files=lambda _:['a'],nanopore_files=lambda _:['b'],read_quast_metric=quast_metric)
    exec(compile(block,str(ROOT/'Hybas.smk'),'exec'),scope)
    with output.open() as f:rows={x['check']:x for x in csv.DictReader(f,delimiter='\t')}
    assert rows['busco_complete_recorded']['status']==expected
    assert rows['analytical_suitability']['status']=='NOT_ASSESSED'
    assert rows['processing_evidence_status']['status'] in {'RECORDED','REVIEW_REQUIRED'}
    assert 'overall_qc_status' not in rows


@pytest.mark.parametrize('changed', ['reads','draft','model','none'])
def test_each_scheduled_medaka_trigger_has_fresh_native_attempt(tmp_path,changed):
    reads,draft,model = [tmp_path/n for n in ('reads.fq','draft.fa','model.tar.gz')]
    fq(reads,['modelA']);draft.write_text('>draft\nAAAA\n');model.write_text('modelA')
    output,receipt=tmp_path/'out/consensus.fa',tmp_path/'out/execution.json';stub=NativeStub(model)
    before=execute(reads,draft,output,receipt,1,requested='modelA',runner=stub,known_models=())
    if changed=='reads':fq(reads,['modelB'])
    elif changed=='draft':draft.write_text('>draft\nCCCC\n')
    elif changed=='model':model.write_text('modelB')
    after=execute(reads,draft,output,receipt,1,requested='modelA',runner=stub,known_models=())
    assert before['attempt']!=after['attempt']
    assert sum(x[0]=='medaka_consensus' for x in stub.calls)==2
    if changed in {'reads','draft'}:assert before['inputs_sha256'][changed]!=after['inputs_sha256'][changed]
    if changed=='model':assert before['model']['sha256']!=after['model']['sha256']


def test_rg_only_metadata_uses_full_model_and_detects_late_mixture(tmp_path):
    reads=tmp_path/'reads.fq'
    reads.write_text('@r RG:Z:run_dna_model@v1.0_barcode01\nACGT\n+\nIIII\n')
    assert model_metadata(reads,('dna_model','dna_model@v1.0'))['basecaller_models']==['dna_model@v1.0']
    with reads.open('a') as f:f.write('@s RG:Z:run_dna_model@v2.0_barcode02\nACGT\n+\nIIII\n')
    assert model_metadata(reads,('dna_model','dna_model@v1.0','dna_model@v2.0'))['status']=='mixed'


def test_resolver_failure_never_runs_consensus_with_native_default(tmp_path):
    reads,draft=tmp_path/'reads.fq',tmp_path/'draft.fa';fq(reads,['modelA']);draft.write_text('>draft\nAAAA\n')
    calls=[]
    def fail(args,**kwargs):
        calls.append(args)
        raise subprocess.CalledProcessError(1,args)
    with pytest.raises(subprocess.CalledProcessError):
        execute(reads,draft,tmp_path/'out/consensus.fa',tmp_path/'out/execution.json',1,runner=fail,known_models=())
    assert len(calls)==1 and calls[0][:3]==['medaka','tools','resolve_model']


def test_v9_absent_platform_confidence_policy_is_preserved(tmp_path):
    c=controller(ROOT/'Hybas.smk')
    left,right=tmp_path/'left',tmp_path/'right'
    left.write_text('100\t100\t0\tG\t10\tExample\n99\t99\t99\tS\t11\tExample one\n');right.write_text('')
    result=c['compute_taxonomic_confidence'](left,right,'11')
    assert result['decision_confidence']=='low'
    assert result['platform_agreement_reason']=='no_inference'
    assert c['choose_platform_cleaning_decision']('aware','11',99,'illumina')[0]=='passthrough'
