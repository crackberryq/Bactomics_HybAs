import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/lib'))
from aware_output_report import decisions, render


def fixture():
    return dict(isolate='sample<script>',limitation='Provisional groups',filtering={'selected_target_taxid':'22','selected_target_name':'Manual choice','kraken_decision_mode':'manual','illumina_kraken_decision':'whitelist_filter','illumina_kraken_decision_reason':'manual target'},parent_recommendation={'name':'Automatic choice','taxid':'11'},statistics=None,groups=[dict(id='primary',name='Primary',category='primary',units={'illumina':3},mode='illumina_only',assembly_eligible=True,reason='supported',manifest=None)])


def test_actual_selection_is_distinct_and_html_is_escaped():
    text=render(fixture())
    assert 'Manual choice · TaxID 22' in text and 'Automatic choice · TaxID 11' in text
    assert 'whitelist_filter' in text and '3 pairs' in text
    assert '<script>' not in text and 'sample&lt;script&gt;' in text
    assert 'not requested' in text and 'does not establish that assembly' in text


def test_statistics_are_readable_and_do_not_change_target():
    d=fixture();d['statistics']=dict(target_taxid='11',target_scope='automatic recommendation',interpretation='Not identity probability',platforms={'illumina':dict(target_units=3,units=4,fraction=.75,wilson_95=[.30,.95])},libraries=[],stability=dict(same_target_fraction=.8,replicates=200,seed=1729,method='stratified bootstrap',selections={'11':160,'22':40}))
    text=render(d)
    assert '75.00%' in text and '30.00% – 95.00%' in text and '80.00%' in text
    assert 'Manual choice · TaxID 22' in text and '<pre>' not in text


def test_decision_sample_and_duplicates_are_checked(tmp_path):
    p=tmp_path/'filtering.tsv';p.write_text('field\tvalue\nisolate\tsample\n')
    assert decisions(p,'sample')['isolate']=='sample'
    with pytest.raises(ValueError):decisions(p,'other')
    p.write_text(p.read_text()+'isolate\tsample\n')
    with pytest.raises(ValueError):decisions(p,'sample')
