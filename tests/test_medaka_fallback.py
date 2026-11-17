"""Documented default selection contracts, not model/biological validation."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
import medaka_stage as stage
from test_medaka_model_review import FAST, HAC, SUP, MAPPING, Native, reads


@pytest.fixture(autouse=True)
def reviewed_runtime(monkeypatch):
    monkeypatch.setattr(stage, 'runtime_identity', lambda: {'medaka_version': '2.2.0',
        'device_policy': {'selection': 'synthetic runtime; no inference'}})
    monkeypatch.setattr(stage, 'installed_defaults', lambda: {'consensus': stage.DOCUMENTED_DEFAULT['model'],
        'options_sha256': stage.DOCUMENTED_DEFAULT['options_sha256']})


def preflight(tmp_path, models, **kwargs):
    native = kwargs.pop('runner', Native(tmp_path))
    path = reads(tmp_path / 'reads.fq', models)
    result = stage.preflight([path], tmp_path / 'receipt.json', model_mapping=MAPPING,
                             runner=native, selection_policy='documented_default', **kwargs)
    return result, native, path


@pytest.mark.parametrize('models,reason', [([None, None], 'basecaller_metadata_absent'),
                                        ([FAST, FAST], 'basecaller_model_unsupported')])
def test_documented_default_is_labelled_and_resolved_once(tmp_path, models, reason, capsys):
    result, native, _ = preflight(tmp_path, models)
    assert result['model_selection_policy'] == 'documented_default'
    assert result['requested_model'] is None and result['model']['requested'] is None
    assert result['input_metadata']['reads'] == len(models)
    assert result['model']['selection'] == 'fallback'
    assert result['model']['selection_reason'] == reason
    assert result['model']['selected_model'] == stage.DOCUMENTED_DEFAULT['model']
    assert result['model']['fallback'] == stage.DOCUMENTED_DEFAULT
    assert result['model']['fallback']['metadata_match'] is False
    assert result['model']['fallback']['scientific_suitability'] == 'not_assessed'
    assert result['model']['fallback']['upstream_source_commit'] == '019a2ab597218a1c239fb7a5315f271a4577ed6b'
    assert len(result['model']['sha256']) == 64
    assert native.calls == [['medaka', 'tools', 'resolve_model', '--model', stage.DOCUMENTED_DEFAULT['model']]]
    assert 'not a metadata match or scientific validation' in capsys.readouterr().out


def test_exact_mapping_and_reviewed_explicit_choice_take_priority(tmp_path, monkeypatch):
    def unexpected_default():
        pytest.fail('A matched or explicit choice must not inspect/substitute the default')
    monkeypatch.setattr(stage, 'installed_defaults', unexpected_default)
    automatic, native, _ = preflight(tmp_path, [HAC])
    assert automatic['model']['selection'] == 'automatic_from_complete_read_metadata'
    assert native.calls[0][-1] == MAPPING[HAC][0]
    explicit = 'custom-research-model'
    selected, native, _ = preflight(tmp_path, [None], requested=explicit,
        review_acknowledged=True, review_reason='Researcher selected this model', reviewed_model=explicit)
    assert selected['model']['selection'] == 'explicit'
    assert selected['model']['requested'] == explicit and native.calls[0][-1] == explicit
    assert 'fallback' not in selected['model']


@pytest.mark.parametrize('models', [[HAC, SUP], [HAC, None], [FAST, None], [FAST, HAC]])
def test_mixed_or_partially_annotated_reads_do_not_select_a_default(tmp_path, models):
    native = Native(tmp_path)
    with pytest.raises(ValueError, match='mixed, incomplete'):
        preflight(tmp_path, models, runner=native)
    assert not native.calls
    assert not (tmp_path / 'receipt.json').exists()


@pytest.mark.parametrize('policy', [None, 'strict'])
def test_omitted_or_strict_policy_keeps_previous_failure_semantics(tmp_path, policy):
    native = Native(tmp_path)
    kwargs = {} if policy is None else {'selection_policy': policy}
    with pytest.raises(ValueError, match='one model identity'):
        stage.preflight([reads(tmp_path / 'reads.fq', [None])], tmp_path / 'receipt.json',
                        model_mapping=MAPPING, runner=native, **kwargs)
    assert not native.calls


@pytest.mark.parametrize('fault', ['wrong_version', 'changed_default', 'changed_source', 'missing_defaults', 'missing_mapping'])
def test_runtime_faults_are_not_treated_as_missing_read_metadata(tmp_path, monkeypatch, fault):
    native = Native(tmp_path)
    if fault == 'wrong_version':
        monkeypatch.setattr(stage, 'runtime_identity', lambda: {'medaka_version': '2.3.0'})
    elif fault == 'changed_default':
        monkeypatch.setattr(stage, 'installed_defaults', lambda: {'consensus': 'changed'})
    elif fault == 'changed_source':
        monkeypatch.setattr(stage, 'installed_defaults', lambda: {'consensus': stage.DOCUMENTED_DEFAULT['model'],
            'options_sha256': 'wrong-source'})
    elif fault == 'missing_defaults':
        def failed(): raise ImportError('broken Medaka installation')
        monkeypatch.setattr(stage, 'installed_defaults', failed)
    with pytest.raises((ValueError, ImportError)):
        stage.preflight([reads(tmp_path / 'reads.fq', [None])], tmp_path / 'receipt.json',
            model_mapping={} if fault == 'missing_mapping' else MAPPING, runner=native,
            selection_policy='documented_default')
    assert not native.calls
    record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
    assert record['status'] == 'failed' and record['model_selection_policy'] == 'documented_default'


def test_corrupt_fastq_after_many_reads_is_not_absent_metadata(tmp_path):
    native = Native(tmp_path)
    path = reads(tmp_path / 'reads.fq', [None] * 101)
    with path.open('a') as handle:
        handle.write('@bad\nACGT\n+\nII\n')
    with pytest.raises(ValueError, match='length mismatch'):
        stage.preflight([path], tmp_path / 'receipt.json', model_mapping=MAPPING, runner=native,
                        selection_policy='documented_default')
    assert not native.calls


@pytest.mark.parametrize('requested,models', [('', [HAC]), ('', [None]), ('chosen-model', [None])])
def test_resolution_or_download_error_never_selects_another_model(tmp_path, requested, models):
    native = Native(tmp_path); native.fail = True
    with pytest.raises(subprocess.CalledProcessError):
        preflight(tmp_path, models, runner=native, requested=requested, review_acknowledged=True,
                  review_reason='Researcher choice', reviewed_model=requested)
    assert len(native.calls) == 1
    record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
    assert record['status'] == 'failed'
    assert record['requested_model'] == (requested or None)
    assert record['resolution_stderr_tail'] == 'original resolver/download failure'
    assert record['model_decision']['selection'] == ('explicit' if requested else
        'fallback' if models == [None] else 'automatic_from_complete_read_metadata')


def test_fallback_executes_medaka_and_preserves_failed_attempt_outputs(tmp_path):
    native = Native(tmp_path)
    result, _, path = preflight(tmp_path, [None], runner=native)
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    def failed_inference(command, **kwargs):
        if command[0] == 'medaka_consensus':
            native.calls.append(command)
            out = Path(command[command.index('-o') + 1]); out.mkdir()
            (out / 'partial.data').write_text('generated intermediate')
            raise subprocess.CalledProcessError(1, command, stderr='device or inference failure')
        return native(command, **kwargs)
    with pytest.raises(subprocess.CalledProcessError):
        stage.execute(path, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 2,
            model_mapping=MAPPING, runner=failed_inference, preflight_receipt=tmp_path / 'receipt.json',
            selection_policy='documented_default')
    assert len(native.calls) == 3  # preflight resolution, execution resolution, one inference; no retries.
    assert not (tmp_path / 'out/consensus.fa').exists()
    partial = list(tmp_path.glob('out/attempts/*/native/partial.data'))
    assert len(partial) == 1 and partial[0].read_text() == 'generated intermediate'
    record = json.loads(next(tmp_path.glob('out/attempts/*/execution.json')).read_text())
    assert record['status'] == 'failed' and record['model']['selection'] == 'fallback'
    assert record['model']['sha256'] == result['model']['sha256']


@pytest.mark.parametrize('field', ['model_selection_policy', 'model_decision'])
def test_execution_binds_selection_policy_and_decision_to_preflight(tmp_path, field):
    result, native, path = preflight(tmp_path, [None])
    receipt = tmp_path / 'receipt.json'
    if field == 'model_selection_policy': result[field] = 'strict'
    else: result[field]['selection_reason'] = 'invented'
    receipt.write_text(json.dumps(result))
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    with pytest.raises(ValueError, match='differs from its preflight'):
        stage.execute(path, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 1,
            model_mapping=MAPPING, runner=native, preflight_receipt=receipt, selection_policy='documented_default')
    assert not any(command[0] == 'medaka_consensus' for command in native.calls)


def test_capability_version_cannot_reuse_the_previous_no_fallback_contract():
    root = Path(__file__).resolve().parents[1]
    marker = json.loads((root / 'app/medaka-model-review.json').read_text())
    assert marker['version'] == stage.MODEL_REVIEW_CONTRACT == 2
    assert marker['missing_policy'] == 'strict'
    assert marker['new_draft_policy'] == 'documented_default'
    for key in ('model', 'medaka_version', 'upstream_source_commit', 'options_sha256'):
        assert marker['automatic_fallback'][key] == stage.DOCUMENTED_DEFAULT[key]


def test_fallback_executes_requested_medaka_stage_and_records_success(tmp_path):
    prior, native, path = preflight(tmp_path, [None])
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    result = stage.execute(path, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 2,
        model_mapping=MAPPING, runner=native, preflight_receipt=tmp_path / 'receipt.json',
        selection_policy='documented_default')
    assert result['status'] == 'completed' and result['model']['selection'] == 'fallback'
    assert result['model_decision'] == prior['model_decision']
    assert result['command'][0] == 'medaka_consensus' and '-m' in result['command']
    assert (tmp_path / 'out/consensus.fa').read_text() == '>synthetic\nACGT\n'
    assert len(list(tmp_path.glob('out/attempts/*/native/consensus.fasta'))) == 1


def test_read_metadata_change_cannot_reclassify_fallback_as_matched_at_execution(tmp_path):
    _, native, _ = preflight(tmp_path, [None])
    filtered = reads(tmp_path / 'filtered.fq', [HAC])
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    # The stub resolves both decisions to the same model bytes. A hash-only
    # comparison would miss the changed decision; the receipt must bind both.
    with pytest.raises(ValueError, match='differs from its preflight'):
        stage.execute(filtered, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 1,
            model_mapping=MAPPING, runner=native, preflight_receipt=tmp_path / 'receipt.json',
            selection_policy='documented_default')
    assert not any(command[0] == 'medaka_consensus' for command in native.calls)
