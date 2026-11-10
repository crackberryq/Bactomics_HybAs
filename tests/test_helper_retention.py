"""File-preservation contracts using synthetic data; no biological tools run."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
import medaka_stage


def test_medaka_failed_publication_scratch_survives_later_attempt(tmp_path, monkeypatch):
    reads, draft, model = [tmp_path / name for name in ('reads.fq', 'draft.fa', 'model')]
    reads.write_text('@fixture model_version_id=model\nACGT\n+\nIIII\n')
    draft.write_text('>draft\nACGT\n')
    model.write_text('synthetic model')
    output, receipt = tmp_path / 'consensus.fa', tmp_path / 'receipt.json'
    old_scratch = tmp_path / 'consensus.fa.tmp'
    old_scratch.write_text('previous incomplete publication')

    def runner(args, **kwargs):
        if args[0] == 'medaka':
            return SimpleNamespace(stdout=str(model) + '\n')
        native = Path(args[args.index('-o') + 1])
        native.mkdir()
        (native / 'consensus.fasta').write_text('>fixture\nACGT\n')
        return SimpleNamespace(returncode=0)

    replace = medaka_stage.os.replace

    def fail_publication(source, destination):
        if Path(destination) == output:
            raise OSError('synthetic publication failure')
        return replace(source, destination)

    with monkeypatch.context() as patch:
        patch.setattr(medaka_stage.os, 'replace', fail_publication)
        with pytest.raises(OSError, match='synthetic publication failure'):
            medaka_stage.execute(reads, draft, output, receipt, 1, runner=runner, known_models=())
    failed = next((tmp_path / 'attempts').iterdir())
    saved = (failed / 'consensus.publish.tmp').read_bytes()
    assert saved == (failed / 'native/consensus.fasta').read_bytes()
    assert json.loads((failed / 'execution.json').read_text())['status'] == 'failed'
    result = medaka_stage.execute(reads, draft, output, receipt, 1, runner=runner, known_models=())
    assert Path(result['attempt']) != failed
    assert (failed / 'consensus.publish.tmp').read_bytes() == saved
    assert old_scratch.read_text() == 'previous incomplete publication'
    assert output.read_bytes() == saved


def test_split_cleaned_library_inputs_survive_preprocessor_failure(tmp_path, monkeypatch):
    import gzip
    import library_inputs
    plan = {'libraries': [dict(library_id='fixture', library_type='paired_end',
                              orientation='fr', read_groups=[{}])]}
    reads = []
    for mate in (1, 2):
        path = tmp_path / f'R{mate}.fq.gz'
        path.write_bytes(gzip.compress(b'@hybasL0G0R1\nACGT\n+\nIIII\n'))
        reads.append(path)

    def fail(*args):
        raise RuntimeError('synthetic preprocessor failure')

    monkeypatch.setattr(library_inputs, 'preprocess_library', fail)
    outdir = tmp_path / 'libraries'
    for _ in range(2):
        with pytest.raises(RuntimeError, match='synthetic preprocessor failure'):
            library_inputs.trim(plan, *reads, outdir, tmp_path / 'out1.gz', tmp_path / 'out2.gz',
                                tmp_path / 'report.html', tmp_path / 'report.json', 1)
    attempts = sorted((outdir / 'cleaned_inputs').iterdir())
    assert len(attempts) == 2
    for attempt in attempts:
        assert list(library_inputs.pairs(attempt / '0_R1.fq.gz', attempt / '0_R2.fq.gz'))
        assert json.loads((attempt / 'failure.json').read_text())['published'] is False


def test_split_cleaned_library_inputs_are_linked_from_completed_metadata(tmp_path):
    import gzip
    import library_inputs
    plan = {'libraries': [dict(library_id='fixture', library_type='paired_end',
                              orientation='fr', preprocessing='pretrimmed', read_groups=[{}])]}
    reads = []
    for mate in (1, 2):
        path = tmp_path / f'R{mate}.fq.gz'
        path.write_bytes(gzip.compress(b'@hybasL0G0R1\nACGT\n+\nIIII\n'))
        reads.append(path)
    library_inputs.trim(plan, *reads, tmp_path / 'libraries', tmp_path / 'out1.gz',
                        tmp_path / 'out2.gz', tmp_path / 'report.html', tmp_path / 'report.json', 1)
    entry = json.loads((tmp_path / 'report.json').read_text())['libraries'][0]
    assert entry['cleaned_pairs'] == entry['trimmed_pairs'] == 1
    assert list(library_inputs.pairs(entry['cleaned_r1'], entry['cleaned_r2']))


from test_aware_preservation import evidence


def test_aware_failed_generated_products_survive_retry(evidence, monkeypatch):
    import aware_preservation
    review, _, output = evidence

    def fail(review, staging, writer_limit):
        (staging / 'generated.fastq').write_text('partial synthetic tool output')
        raise RuntimeError('synthetic preservation failure')

    monkeypatch.setattr(aware_preservation, '_preserve', fail)
    for _ in range(2):
        with pytest.raises(RuntimeError, match='synthetic preservation failure'):
            aware_preservation.preserve(review, output)
    attempts = sorted(output.parent.glob(output.name + '.attempt-*'))
    assert len(attempts) == 2 and not output.exists()
    for attempt in attempts:
        assert (attempt / 'generated.fastq').read_text() == 'partial synthetic tool output'
        assert json.loads((attempt / 'failure.json').read_text())['published'] is False
    assert not output.with_name(output.name + '.aware-plus.lock').exists()


def test_aware_bundle_failed_primary_read_products_are_retained(evidence, monkeypatch):
    import aware_bundle
    import aware_output_report
    review, paths, output = evidence
    decision = output.parent / 'decision.tsv'
    decision.write_text('field\tvalue\nisolate\tisolate1\n')

    def fail(document):
        raise RuntimeError('synthetic report failure')

    monkeypatch.setattr(aware_output_report, 'render', fail)
    with pytest.raises(RuntimeError, match='synthetic report failure'):
        aware_bundle.build(review, output, {'illumina': [paths['r1'], paths['r2']]}, decision)
    attempt = next(output.parent.glob(output.name + '.attempt-*'))
    assert list((attempt / 'primary').rglob('*.fastq.gz'))
    assert json.loads((attempt / 'failure.json').read_text())['published'] is False
    assert not output.exists()
    with pytest.raises(ValueError, match='unpublished failed attempt'):
        aware_bundle.verify_bundle(attempt / 'manifest.json')


def test_aware_completed_bundle_inventories_retained_bootstrap_reports(evidence):
    from aware_bundle import build, verify_bundle
    review, paths, output = evidence
    receipt = json.loads(review.read_text())
    receipt['recommendation'] = {'taxid': '11'}
    review.write_text(json.dumps(receipt))
    decision = output.parent / 'decision.tsv'
    decision.write_text('field\tvalue\nisolate\tisolate1\n')
    build(review, output, {'illumina': [paths['r1'], paths['r2']]}, decision,
          statistics=True, selector=lambda reports: {'selected_taxid': '11'})
    document = verify_bundle(output / 'manifest.json')
    reports = [row for row in document['files'] if row['path'].startswith('bootstrap_reports/')]
    assert len(reports) == 600
    assert all((output / row['path']).is_file() for row in reports)
    assert not (output / 'failure.json').exists()


def test_bootstrap_reports_are_retained_per_replicate_without_changing_results(evidence):
    from aware_statistics import analyse
    review, _, output = evidence
    receipt = json.loads(review.read_text())
    receipt['recommendation'] = {'taxid': '11'}
    snapshots = {}

    def selector(reports):
        for path in reports.values():
            snapshots[path] = path.read_bytes()
        return {'selected_taxid': '11'}

    first = analyse(receipt, selector, replicates=4, report_dir=output)
    second = analyse(receipt, selector, replicates=4, report_dir=output)
    assert first == second
    assert first['stability']['selections'] == {'11': 4}
    assert len(snapshots) == 24
    assert all(path.read_bytes() == contents for path, contents in snapshots.items())


def test_bootstrap_failed_selector_keeps_generated_reports(evidence):
    from aware_statistics import analyse
    review, _, output = evidence
    receipt = json.loads(review.read_text())
    receipt['recommendation'] = {'taxid': '11'}
    observed = []

    def fail(reports):
        observed.extend(reports.values())
        raise RuntimeError('synthetic selector failure')

    with pytest.raises(RuntimeError, match='synthetic selector failure'):
        analyse(receipt, fail, replicates=4, report_dir=output)
    assert len(observed) == 3 and all(path.is_file() for path in observed)
    assert any(path.stat().st_size for path in observed)


def test_computational_report_refresh_retains_previous_document(tmp_path):
    from computational import render_report
    report = tmp_path / "computational_report.html"
    report.write_text("previous computational document")
    render_report(tmp_path)
    saved = list((tmp_path / "report_history").glob("*.html"))
    assert len(saved) == 1
    assert saved[0].read_text() == "previous computational document"
    current = report.read_bytes()
    render_report(tmp_path)
    assert current in [p.read_bytes() for p in (tmp_path / "report_history").glob("*.html")]
