"""Software-only preservation fixtures, not strain-resolution validation."""
from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/lib'))
from aware_preservation import FILES, preserve
from identification_review import identity, SCHEMA as REVIEW_SCHEMA

NODES = [('1','R','root',0,None),('2','D','Bacteria',2,'1'),
         ('10','G','Example',4,'2'),('11','S','Example alpha',6,'10'),
         ('111','S1','Example alpha strain one',8,'11'),
         ('12','S','Example beta',6,'10')]


def report(path, calls):
    direct=Counter(r[1] for r in calls)
    total=sum(direct.values());clades=direct.copy()
    for taxid,rank,name,depth,parent in reversed(NODES):
        if parent:clades[parent]+=clades[taxid]
    rows=[]
    if direct['0']:rows.append(f'{100*direct["0"]/total}\t{direct["0"]}\t{direct["0"]}\tU\t0\tunclassified')
    for taxid,rank,name,depth,parent in NODES:
        rows.append(f'{100*clades[taxid]/total if total else 0}\t{clades[taxid]}\t{direct[taxid]}\t{rank}\t{taxid}\t'+(' '*depth)+name)
    path.write_text('\n'.join(rows)+'\n')


def write_reads(path, names, mate=None):
    records=[]
    for name in names:
        label=name+('/'+str(mate) if mate else '')
        records.append('@'+label+' fixture\nACGTACGT\n+\nIIIIIIII\n')
    path.write_bytes(gzip.compress(''.join(records).encode(),mtime=0))


def receipt(path, paths, illumina=True, ont=True):
    path.write_text(json.dumps({'schema':REVIEW_SCHEMA,'isolate':'isolate1',
        'context':'artificial-review-context','platforms':{'illumina':illumina,'ont':ont},
        'files':[identity(p) for p in paths.values()]}))


@pytest.fixture
def evidence(tmp_path):
    folder=tmp_path/'input evidence';folder.mkdir()
    paths={k:folder/v for k,v in FILES.items()}
    short=[('pair1','11','11:3 |:| 11:3'),('pair2','12','12:3 |:| 12:3'),
           ('pair3','10','10:3 |:| 10:3'),('pair4','0','0:3 |:| 0:3')]
    long=[('long1','11','11:3 10:1'),('long2','12','12:4'),
          ('long3','10','10:4'),('long4','11','11:2 12:2'),
          ('long5','11','11:2 999:2'),('long6','0','0:4'),('long7','111','111:4')]
    for platform,calls in [('illumina',short),('ont',long)]:
        report(paths[platform+'_report'],calls)
        paths[platform+'_calls'].write_text(''.join(
            f'{"U" if taxid=="0" else "C"}\t{name}\t{taxid}\t{"8|8" if platform=="illumina" else "8"}\t{hits}\n'
            for name,taxid,hits in reversed(calls)))
    write_reads(paths['r1'],[r[0] for r in short],1)
    write_reads(paths['r2'],[r[0] for r in short],2)
    write_reads(paths['ont'],[r[0] for r in long])
    p=folder/'review.json';receipt(p,paths)
    return p,paths,tmp_path/'preserved output'


def test_preserves_all_units_with_pairs_and_shared_pool(evidence):
    review,paths,output=evidence
    before={p:identity(p) for p in paths.values()}
    result=preserve(review,output,writer_limit=1)
    assert result['input_units']==result['output_units']=={'illumina':4,'ont':7}
    bins={b['taxid'] or b['category']:b for b in result['bins']}
    assert bins['11']['units']=={'illumina':1,'ont':2}
    assert bins['12']['units']=={'illumina':1,'ont':1}
    assert bins['ambiguous']['units']=={'illumina':1,'ont':3}
    assert bins['unclassified']['units']=={'illumina':1,'ont':1}
    assert all(b['reconstruction_status']=='not_available' for b in result['bins'])
    with gzip.open(output/'assignments.tsv.gz','rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==11 and len({(r['platform'],r['read_id']) for r in rows})==11
    shared=next(r for r in rows if r['read_id']=='long3')
    assert shared['bucket']=='_ambiguous' and shared['candidate_species_taxids']=='11,12'
    conflict=next(r for r in rows if r['read_id']=='long4')
    assert conflict['reason']=='conflicting_species_hits'
    unknown=next(r for r in rows if r['read_id']=='long5')
    assert unknown['unresolved_taxonomy']=='True'
    group=output/bins['11']['directory']
    assert group.name=='isolate1_example_alpha_taxid_11'
    assert gzip.decompress((group/'illumina_L0/G0_r1.fastq.gz').read_bytes()).startswith(b'@pair1/1 ')
    assert gzip.decompress((group/'illumina_L0/G0_r2.fastq.gz').read_bytes()).startswith(b'@pair1/2 ')
    assert len(gzip.decompress((group/'ont_L0/G0_reads.fastq.gz').read_bytes()).splitlines())==8
    assert result['metrics']['wall_seconds']>0 and result['metrics']['peak_rss_kib']>0
    assert result['review_context']=='artificial-review-context'
    for p,original in before.items():assert identity(p)==original
    for item in result['outputs']:
        actual=identity(output/item['path']);assert actual['sha256']==item['sha256']
    assert 'not reconstructed isolates' in (output/'report.html').read_text()
    assert not list(output.glob('*.sqlite*'))


def test_ont_only_keeps_absent_platform_empty(evidence):
    review,paths,output=evidence
    for key in ('r1','r2','illumina_calls','illumina_report'):paths[key].write_bytes(b'')
    receipt(review,paths,illumina=False)
    result=preserve(review,output)
    assert result['input_units']=={'ont':7}
    assert not list(output.rglob('illumina_*.fastq.gz'))


def test_short_reads_do_not_require_long_reads(evidence):
    review,paths,output=evidence
    for key in ('ont','ont_calls','ont_report'):paths[key].write_bytes(b'')
    receipt(review,paths,ont=False)
    result=preserve(review,output)
    assert result['input_units']==result['output_units']=={'illumina':4}


def test_stale_evidence_is_rejected(evidence):
    review,paths,output=evidence
    paths['ont_calls'].write_text(paths['ont_calls'].read_text()+'\n')
    with pytest.raises(ValueError,match='Reviewed input changed'):preserve(review,output)
    assert not output.exists()


@pytest.mark.parametrize('change, message', [
    ('duplicate_call','Duplicate Kraken'),('duplicate_read','duplicate FASTQ'),
    ('missing_read','absent from FASTQ'),('missing_call','counts do not match'),
    ('wrong_length','lengths differ'),('wrong_mate','mate identifiers'),
    ('extra_mate','mate counts'),('malformed_fastq','Malformed'),
    ('bad_hit','hit list'),('bad_pair_separator','separator'),
    ('report_disagreement','disagree'),('wrong_report_count','counts do not match')])
def test_invalid_inputs_do_not_publish_partial_results(evidence,change,message):
    review,paths,output=evidence
    if change=='duplicate_call':
        paths['ont_calls'].write_text(paths['ont_calls'].read_text()*2)
    elif change=='duplicate_read':write_reads(paths['ont'],['long1','long1','long2','long3','long4','long5','long6','long7'])
    elif change=='missing_read':write_reads(paths['ont'],['long1','long2','long3','long4','long5','long6'])
    elif change=='missing_call':paths['ont_calls'].write_text('\n'.join(paths['ont_calls'].read_text().splitlines()[1:])+'\n')
    elif change=='wrong_length':paths['ont_calls'].write_text(paths['ont_calls'].read_text().replace('\t8\t','\t9\t'))
    elif change=='wrong_mate':write_reads(paths['r2'],['other','pair2','pair3','pair4'],2)
    elif change=='extra_mate':write_reads(paths['r2'],['pair1','pair2','pair3'],2)
    elif change=='malformed_fastq':paths['ont'].write_bytes(gzip.compress(b'@long1\nACGT\n+\nI\n'))
    elif change=='bad_hit':paths['ont_calls'].write_text(paths['ont_calls'].read_text().replace('0:4','broken:4'))
    elif change=='bad_pair_separator':paths['illumina_calls'].write_text(paths['illumina_calls'].read_text().replace('|:|','|:| |:|'))
    elif change=='report_disagreement':paths['ont_report'].write_text(paths['ont_report'].read_text().replace('Example alpha\n','Different alpha\n'))
    elif change=='wrong_report_count':paths['ont_report'].write_text(paths['ont_report'].read_text().replace('\t1\tS1\t111\t','\t2\tS1\t111\t'))
    receipt(review,paths)
    untouched=output.parent/'keep.txt';untouched.write_text('preserve me')
    with pytest.raises(ValueError,match=message):preserve(review,output)
    assert not output.exists() and untouched.read_text()=='preserve me'
    assert not list(output.parent.glob('.preserved output-*'))
    assert not output.with_name(output.name+'.aware-plus.lock').exists()


def test_existing_output_and_concurrent_writer_are_retained(evidence):
    review,paths,output=evidence
    output.mkdir();marker=output/'keep.txt';marker.write_text('saved results')
    with pytest.raises(ValueError,match='already exists'):preserve(review,output)
    assert marker.read_text()=='saved results'
    other=output.parent/'other';lock=other.with_name(other.name+'.aware-plus.lock');lock.write_text('active fixture')
    with pytest.raises(FileExistsError):preserve(review,other)
    assert lock.read_text()=='active fixture'
