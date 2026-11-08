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
