"""Regression contract against the exact active v8.8 controller (no biological inputs)."""
import importlib.util
import os
from pathlib import Path
import random
import re
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflow_support import manual_taxid, normalize_taxid, pair_illumina_files, polish_plan, validate_isolate_names, busco_db_ready, quast_metric

DEFAULTS = dict(KRAKEN_AUTO_PREFER_SPECIES_MIN_PCT=90., KRAKEN_AUTO_PASSTHROUGH_MIN_PCT=95.,
    KRAKEN_AUTO_USE_BOTH_PLATFORMS=True, KRAKEN_AWARE_DELTA_MIN_PCT=20.,
    KRAKEN_AWARE_MODERATE_DELTA_MIN_PCT=10., KRAKEN_AWARE_HIGH_TARGET_MIN_PCT=95.,
    KRAKEN_AWARE_MODERATE_TARGET_MIN_PCT=80.)


def controller(path):
    text = path.read_text()
    names = {'normalize_taxid', 'safe_float', 'combine_platform_values', 'choose_platform_cleaning_decision',
        'parse_kraken_report_rows', 'parse_kraken_report_target_fraction', 'is_meaningful_kraken_candidate',
        'infer_target_from_kraken_report', 'choose_selected_kraken_target', 'choose_kraken_decision',
        'top_two_meaningful_candidates', 'compute_platform_agreement', 'compute_taxonomic_confidence'}
    scope = dict(os=os, normalize_taxid=normalize_taxid, **DEFAULTS)
    for block in re.findall(r'^def [\s\S]*?(?=^\S|\Z)', text, re.M):
        if block.split('(', 1)[0].removeprefix('def ') in names:
            exec(compile(block, str(path), 'exec'), scope)
    return scope


def evaluate(c, mode, target, left, right):
    selected = c['choose_selected_kraken_target'](mode, target, left, right)
    return (selected, c['choose_kraken_decision'](mode, selected['selected_taxid'], selected['selected_pct_for_decision']),
        c['choose_platform_cleaning_decision'](mode, selected['selected_taxid'], selected['illumina_selected_pct'], 'illumina'),
        c['choose_platform_cleaning_decision'](mode, selected['selected_taxid'], selected['ont_selected_pct'], 'ont'),
        c['compute_taxonomic_confidence'](left, right, selected['selected_taxid']))


def test_5000_default_policy_cases(tmp_path):
    old = controller(ROOT / 'tests/fixtures/aware_v8_8.py')
    new = controller(ROOT / 'Hybas.smk')
    rng = random.Random(9000)
    left, right = tmp_path/'left.tsv', tmp_path/'right.tsv'
    for case in range(5000):
        for path in (left, right):
            rows = []
            if rng.random() > .15:
                for taxid, rank, name in [('11','G','Example'),('12','S','Example one'),('21','G','Other'),('22','S','Other two')]:
                    percent = rng.choice([0, 5, 79.99, 80, 89.99, 90, 94.99, 95, 100, rng.uniform(0,100)])
                    rows.append(f'{percent:.4f}\t100\t50\t{rank}\t{taxid}\t{name}\n')
            path.write_text(''.join(rows))
        mode = rng.choice(['off','manual','auto','aware'])
        target = rng.choice(['','11','12','21','22','999'])
        both = rng.choice([True, False])
        old['KRAKEN_AUTO_USE_BOTH_PLATFORMS'] = new['KRAKEN_AUTO_USE_BOTH_PLATFORMS'] = both
        assert evaluate(old, mode, target, left, right) == evaluate(new, mode, target, left, right), case


def test_confidence_config_is_effective_and_does_not_change_platform_action(tmp_path):
    c = controller(ROOT/'Hybas.smk')
    left, right = tmp_path/'l', tmp_path/'r'
    left.write_text('90\t90\t90\tG\t11\tExample\n')
    right.write_text(left.read_text())
    assert c['compute_taxonomic_confidence'](left,right,'11')['decision_confidence'] == 'moderate'
    before = c['choose_platform_cleaning_decision']('aware','11',90,'illumina')
    c['KRAKEN_AWARE_HIGH_TARGET_MIN_PCT'] = 85
    assert c['compute_taxonomic_confidence'](left,right,'11')['decision_confidence'] == 'high'
    assert c['choose_platform_cleaning_decision']('aware','11',90,'illumina') == before


@pytest.mark.parametrize('value', [None,'','  ',123,'123','123.0',123.0])
def test_manual_taxid_valid(value):
    assert manual_taxid(value) == ('' if value is None or not str(value).strip() else '123')


@pytest.mark.parametrize('value', ['None','null','NaN','abc',-1,True,False,0,'0','1.5','12; echo x'])
def test_manual_taxid_rejects_invalid(value):
    with pytest.raises(ValueError): manual_taxid(value)


@pytest.mark.parametrize('names', [['a','a'],['../a'],['/a'],['a b'],['a;ls'],['']])
def test_unsafe_or_duplicate_names(names):
    with pytest.raises(ValueError): validate_isolate_names(names)


def test_pairing_rejects_mismatched_lane_and_orphans():
    with pytest.raises(ValueError): pair_illumina_files(['x_L1_R1.fastq.gz'], ['x_L2_R2.fastq.gz'])
    with pytest.raises(ValueError): pair_illumina_files(['x_R1.fastq.gz'], [])
    assert pair_illumina_files(['x_L2_R1.fq','x_L1_R1.fastq.gz'], ['x_L1_R2.fq.gz','x_L2_R2.fastq']) == (['x_L1_R1.fastq.gz','x_L2_R1.fq'],['x_L1_R2.fq.gz','x_L2_R2.fastq'])


@pytest.mark.parametrize('mode', ['hybrid','illumina_only','ont_only'])
@pytest.mark.parametrize('requested', ['none','racon','medaka','tripolish'])
@pytest.mark.parametrize('enabled', [False,True])
@pytest.mark.parametrize('rounds', [0,2])
def test_polish_routes(mode,requested,enabled,rounds):
    plan = polish_plan(mode,requested,rounds,enabled)
    if mode == 'illumina_only':
        assert plan == (['polypolish'] if enabled and requested=='tripolish' else [])
    elif requested == 'none': assert plan == []
    else:
        assert plan[:rounds] == [f'racon{i}' for i in range(1,rounds+1)]
        assert ('medaka' in plan) == (requested in {'medaka','tripolish'})
        assert ('polypolish' in plan) == (mode=='hybrid' and enabled and requested=='tripolish')
    if mode=='hybrid' and requested=='tripolish' and not enabled: assert plan[-1]=='medaka'


def test_incomplete_busco_directory_is_not_ready(tmp_path):
    (tmp_path/'.complete').touch()
    assert not busco_db_ready(tmp_path)
    (tmp_path/'hmms').mkdir()
    for path in ['dataset.cfg','scores_cutoff','lengths_cutoff','hmms/test.hmm']: (tmp_path/path).write_text('fixture\n')
    assert busco_db_ready(tmp_path)


def test_quast_native_and_transposed_tables(tmp_path):
    path=tmp_path/'report.tsv'
    path.write_text('Assembly\tassembly\n# contigs\t3\nN50\t1000\n')
    assert quast_metric(path,'N50')=='1000'
    path.write_text('Assembly\t# contigs\tN50\nassembly\t3\t1000\n')
    assert quast_metric(path,'N50')=='1000'
