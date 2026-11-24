"""Saved-output presentation, without biological analysis or native tools."""
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/lib'))
from hybas_report import collect, render_report, normalize, write_report


def fixture(root):
    (root/'metadata').mkdir(parents=True)
    summary=dict(isolate='sample',workflow_version='HybAs-9.0.0-AWARE',run_mode='hybrid',assembly_bp='4529685',
        assembly_contig_count_fasta='1',quast_n50='4529685',quast_gc_percent='38.5',run_busco=True,run_prokka=True,
        busco_lineage='example_odb10',busco_complete_percent='98.4',busco_single_copy_percent='97.4',busco_duplicated_percent='1',busco_fragmented_percent='0.6',busco_missing_percent='1',
        illumina_raw_r1_reads='1000',illumina_clean_r1_reads='900',illumina_trim_r1_reads='850',
        illumina_raw_r2_reads='1000',illumina_clean_r2_reads='900',illumina_trim_r2_reads='850',
        ont_raw_reads='200',ont_clean_reads='190',ont_filt_reads='150',
        kraken_selected_target_name='Example target',kraken_selected_target_taxid='12',kraken_selected_target_source='illumina',
        kraken_illumina_selected_target_percent='99',kraken_ont_selected_target_percent='40',kraken_decision_confidence='low',
        executed_stages=['racon1','racon2','medaka'],skipped_stages={'polypolish':'disabled'},final_source_stage='medaka',
        requested_polish_mode='tripolish',tool_unicycler_version='0.5.1',estimated_ont_coverage_x='35.7')
    (root/'metadata/run_summary.json').write_text(json.dumps(summary))
    (root/'metadata/kraken_decision.tsv').write_text('metric\tvalue\nillumina_kraken_decision\tpass_through\nont_kraken_decision\textract\n')
    (root/'annotation').mkdir();(root/'annotation/sample.gff').write_text('##gff-version 3\nc\tfixture\tCDS\t1\t3\t.\t+\t0\tID=p\n##FASTA\n>c\nAAA\n')
    return summary


def test_collect_roundtrip_preserves_source_and_separates_decisions(tmp_path):
    s=fixture(tmp_path);before=(tmp_path/'metadata/run_summary.json').read_bytes()
    write_report(tmp_path,tmp_path/'reports/hybas')
    data=json.loads((tmp_path/'reports/hybas/results.json').read_text());html=(tmp_path/'reports/hybas/report.html').read_text()
    assert data['summary']==s and data['annotation']['counts']=={'CDS':1}
    assert 'pass_through' in html and 'extract' in html and 'low' in html
    assert '4,529,685' in html and 'Not recorded' in html
    assert data['sources_sha256']['metadata/run_summary.json']
    assert (tmp_path/'metadata/run_summary.json').read_bytes()==before
    assert (tmp_path/'reports/hybas/methods.txt').is_file()


def test_legacy_missing_and_zero_are_distinct():
    s=dict(isolate='legacy',workflow_version='old',assembly_bp=0,run_busco=False)
    text=render_report(s)
    assert 'legacy run summary' in text and 'Disabled in the recorded configuration' in text
    assert '<div class="stat-value">0</div>' in text
    assert '<div class="stat-value">Not recorded</div>' in text
    with pytest.raises(ValueError): normalize({})


def test_untrusted_labels_are_text_not_html():
    text=render_report(dict(isolate='<script>alert(1)</script>',workflow_version='old',kraken_selected_target_name='<img src=x onerror=alert(1)>'))
    assert '<script>alert(1)</script>' not in text and '<img src=x' not in text
    assert '&lt;script&gt;' in text


def test_pacbio_report_retains_platform_and_assembler():
    text=render_report(dict(isolate='sample',workflow_version='10',run_mode='pacbio_only',
        long_read_platform='pacbio',assembler='flye',pacbio={'raw':{'reads':123,'bases':456},'cleaned':{'reads':100,'bases':400}},pacbio_selected_target_percent=88))
    assert 'PacBio' in text and 'flye' in text and '123 reads' in text
    assert 'Oxford Nanopore' not in text and 'After Filtlong' not in text and 'unicycler' not in text


def test_processing_checks_are_not_analytical_suitability():
    data=dict(schema_version='hybas-report-1',summary=dict(isolate='sample',workflow_version='10',run_busco=True,busco_complete_percent=0,executed_stages=['medaka'],medaka_model_requested='reviewed',medaka_model_effective='effective.tar.gz',medaka_model_selection='explicit',medaka_input_metadata_status='mixed'),
              tables={'qc_assessment':[{'check':'overall_qc_status','status':'PASS','detail':'legacy bookkeeping'}]})
    text=render_report(data)
    assert 'Not assessed against validated use-specific acceptance criteria' in text
    assert 'legacy overall_qc_status is also a bookkeeping check' in text
    assert 'Medaka model provenance' in text and 'effective.tar.gz' in text and 'mixed' in text


def test_unmatched_model_review_is_visible_and_escaped():
    text = render_report(dict(isolate='sample', workflow_version='10.0.0a5', executed_stages=['medaka'],
        medaka_model_effective='explicit-model', medaka_model_review_status='unsupported_basecaller',
        medaka_model_review_acknowledged=True, medaka_model_reviewed_model='explicit-model',
        medaka_model_review_reason='<script>research comparison</script>'))
    assert 'unsupported_basecaller' in text and 'Model bound to review' in text
    assert '&lt;script&gt;research comparison&lt;/script&gt;' in text
    assert '<script>research comparison</script>' not in text
    assert 'does not validate model suitability' in text


def test_fallback_is_visible_in_report_and_methods_without_claiming_a_match():
    from hybas_report import methods_text
    summary = dict(isolate='sample', workflow_version='development', executed_stages=['medaka'],
        medaka_model_requested='', medaka_model_selection='fallback',
        medaka_model_selection_policy='documented_default',
        medaka_model_selection_reason='basecaller_metadata_absent',
        medaka_model_fallback_default='r1041_e82_400bps_sup_v5.2.0',
        medaka_model_fallback_medaka_version='2.2.0',
        medaka_model_fallback_source_commit='019a2ab597218a1c239fb7a5315f271a4577ed6b',
        medaka_model_effective='resolved.tar.gz', medaka_model_sha256='a' * 64)
    for text in (render_report(summary), methods_text({'summary': summary})):
        assert 'basecaller_metadata_absent' in text
        assert 'r1041_e82_400bps_sup_v5.2.0' in text
        assert 'not a match to the read metadata' in text
        assert 'does not establish scientific suitability' in text


def test_historical_off_mode_overrides_legacy_confidence_in_presentation():
    s = dict(isolate='sample', workflow_version='10', kraken_decision_mode='off',
             kraken_awareness_state='fail_no_target', kraken_decision_confidence='low')
    page = render_report(s)
    visible = page.split('<script type="application/json"', 1)[0]
    assert 'Kraken skipped by configuration' in visible
    assert 'fail_no_target' not in visible
    assert s['kraken_awareness_state'] == 'fail_no_target'
