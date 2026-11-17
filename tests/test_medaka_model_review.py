"""Model-decision contracts with synthetic reads, never biological validation."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
import medaka_stage as stage

# Identifiers/mappings from nanoporetech/medaka v2.2.0, medaka/options.py.
HAC = 'dna_r10.4.1_e8.2_400bps_hac@v4.2.0'
SUP = 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0'
FAST = 'dna_r10.4.1_e8.2_400bps_fast@v4.2.0'
MAPPING = {HAC: ('r1041_e82_400bps_hac_v4.2.0', None), SUP: ('r1041_e82_400bps_sup_v4.2.0', None)}


def reads(path, models, key='basecall_model_version_id'):
    path.write_text(''.join(f'@r{i}' + (f' {key}={model}' if model else '') + '\nACGT\n+\nIIII\n'
                            for i, model in enumerate(models)))
    return path


class Native:
    def __init__(self, folder):
        self.calls = []
        self.model = folder / 'synthetic-model.tar.gz'
        self.model.write_text('synthetic model; no inference')
        self.fail = False

    def __call__(self, command, **kwargs):
        self.calls.append(command)
        if self.fail:
            raise subprocess.CalledProcessError(1, command, stderr='original resolver/download failure')
        if command[0] == 'medaka':
            return SimpleNamespace(stdout=str(self.model) + '\n')
        output = Path(command[command.index('-o') + 1]); output.mkdir()
        (output / 'consensus.fasta').write_text('>synthetic\nACGT\n')
        return SimpleNamespace(returncode=0)


@pytest.mark.parametrize('model', [HAC, SUP])
def test_mapping_matched_explicit_needs_no_acknowledgement(tmp_path, model):
    native = Native(tmp_path)
    result = stage.preflight([reads(tmp_path / 'reads.fq', [model])], tmp_path / 'receipt.json',
                             requested=MAPPING[model][0], model_mapping=MAPPING, runner=native)
    assert result['model_review']['status'] == 'matched_to_installed_mapping'
    assert not result['model_review']['acknowledgement_required']
    assert result['model_review']['scientific_suitability'] == 'not_assessed'


@pytest.mark.parametrize('models,requested,status', [
    ([FAST], MAPPING[SUP][0], 'unsupported_basecaller'),
    ([HAC], MAPPING[SUP][0], 'mismatched_to_installed_mapping'),
    ([HAC], '/custom/model.tar.gz', 'unverified_explicit_model'),
    ([HAC] * 100 + [SUP], MAPPING[HAC][0], 'metadata_review_required'),
    ([None], MAPPING[SUP][0], 'metadata_review_required'),
])
def test_unmatched_requires_acknowledgement_and_reason_then_stays_labelled(tmp_path, models, requested, status, capsys):
    native = Native(tmp_path); path = reads(tmp_path / 'reads.fq', models)
    receipt = tmp_path / 'receipt.json'
    for ack, reason in [(False, ''), (True, ''), (False, 'Research decision')]:
        with pytest.raises(ValueError, match='requires review'):
            stage.preflight([path], receipt, requested=requested, model_mapping=MAPPING,
                            review_acknowledged=ack, review_reason=reason, reviewed_model=requested, runner=native)
        assert not native.calls and not receipt.exists()
    result = stage.preflight([path], receipt, requested=requested, model_mapping=MAPPING,
                             review_acknowledged=True, review_reason='Reviewed research configuration, not a supported match', reviewed_model=requested, runner=native)
    assert result['model_review']['status'] == status
    assert result['model_review']['review_reason'].startswith('Reviewed research')
    assert 'scientific suitability is not established' in capsys.readouterr().out
    assert result['input_metadata']['reads'] == len(models)
    assert result['input_metadata']['inspection'] == 'all_reads'


def test_failure_receipt_keeps_request_metadata_and_original_error(tmp_path):
    native = Native(tmp_path); native.fail = True
    path = reads(tmp_path / 'reads.fq', [HAC], 'basecaller_model_version_id')
    with pytest.raises(subprocess.CalledProcessError):
        stage.preflight([path], tmp_path / 'receipt.json', model_mapping=MAPPING, runner=native)
    record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
    assert record['requested_model'] is None
    assert record['input_metadata']['basecaller_models'] == [HAC]
    assert record['model_review']['status'] == 'matched_to_installed_mapping'
    assert record['resolution_command'][:3] == ['medaka', 'tools', 'resolve_model']
    assert record['resolution_stderr_tail'] == 'original resolver/download failure'
    assert record['status'] == 'failed'
    assert len(native.calls) == 1


def test_preflight_inspects_every_file_and_execution_rechecks_filtered_reads(tmp_path):
    native = Native(tmp_path)
    first = reads(tmp_path / 'first.fq', [HAC] * 100)
    second = reads(tmp_path / 'second.fq', [SUP])
    receipt = tmp_path / 'receipt.json'
    with pytest.raises(ValueError, match='one model identity'):
        stage.preflight([first, second], receipt, model_mapping=MAPPING, runner=native)
    record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
    assert record['input_metadata']['reads'] == 101
    assert len(record['input_metadata']['files']) == 2 and not native.calls
    stage.preflight([first], receipt, model_mapping=MAPPING, runner=native)
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    filtered = reads(tmp_path / 'filtered.fq', [HAC] * 100 + [SUP])
    with pytest.raises(ValueError, match='one model identity'):
        stage.execute(filtered, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 1,
                      model_mapping=MAPPING, runner=native, preflight_receipt=receipt)
    assert not any(command[0] == 'medaka_consensus' for command in native.calls)


def test_model_replacement_after_preflight_never_runs_consensus(tmp_path):
    native = Native(tmp_path); path = reads(tmp_path / 'reads.fq', [HAC])
    receipt = tmp_path / 'receipt.json'
    stage.preflight([path], receipt, model_mapping=MAPPING, runner=native)
    native.model.write_text('changed synthetic model')
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    with pytest.raises(ValueError, match='differs from its preflight'):
        stage.execute(path, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 1,
                      model_mapping=MAPPING, runner=native, preflight_receipt=receipt)
    assert not any(command[0] == 'medaka_consensus' for command in native.calls)


def test_resolution_command_is_durable_before_native_process(tmp_path):
    native = Native(tmp_path); path = reads(tmp_path / 'reads.fq', [HAC])
    def inspect(command, **kwargs):
        record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
        assert record['resolution_command'] == command
        assert record['input_metadata']['inspection'] == 'all_reads'
        return native(command, **kwargs)
    stage.preflight([path], tmp_path / 'receipt.json', model_mapping=MAPPING, runner=inspect)


@pytest.mark.parametrize('basecaller', [FAST, HAC.replace('@v4.2.0', '@v99.0.0'), HAC + '.1'])
def test_unknown_rg_version_cannot_fall_back_to_unversioned_model(tmp_path, basecaller):
    mapping = {**MAPPING, HAC.split('@')[0]: ('r1041_e82_400bps_hac_g632', None)}
    path = tmp_path / 'reads.fq'
    path.write_text(f'@r RG:Z:run_{basecaller}_barcode01\nACGT\n+\nIIII\n')
    native = Native(tmp_path)
    with pytest.raises(ValueError, match='unsupported_basecaller'):
        stage.preflight([path], tmp_path / 'receipt.json', model_mapping=mapping, runner=native)
    assert not native.calls
    record = json.loads(next(tmp_path.glob('preflight_attempts/*/execution.json')).read_text())
    assert record['input_metadata']['basecaller_models'] == [basecaller]


def test_automatic_resolution_uses_exact_all_read_mapping(tmp_path):
    native = Native(tmp_path)
    stage.preflight([reads(tmp_path / 'reads.fq', [HAC])], tmp_path / 'receipt.json',
                    model_mapping=MAPPING, runner=native)
    assert native.calls == [['medaka', 'tools', 'resolve_model', '--model', MAPPING[HAC][0]]]


@pytest.mark.parametrize('field,value', [('schema', 'unrelated'), ('review_reason', 'different decision'),
                                        ('review_acknowledged', True)])
def test_execution_rejects_unrelated_preflight_decision(tmp_path, field, value):
    native = Native(tmp_path); path = reads(tmp_path / 'reads.fq', [HAC])
    receipt = tmp_path / 'receipt.json'
    prior = stage.preflight([path], receipt, model_mapping=MAPPING, runner=native)
    if field == 'schema': prior[field] = value
    else: prior['model_review'][field] = value
    receipt.write_text(json.dumps(prior))
    draft = tmp_path / 'draft.fa'; draft.write_text('>draft\nACGT\n')
    with pytest.raises(ValueError, match='differs from its preflight'):
        stage.execute(path, draft, tmp_path / 'out/consensus.fa', tmp_path / 'out/receipt.json', 1,
                      model_mapping=MAPPING, runner=native, preflight_receipt=receipt)
    assert not any(command[0] == 'medaka_consensus' for command in native.calls)


@pytest.mark.parametrize('binding', ['', 'another-model'])
def test_unmatched_acknowledgement_is_bound_to_requested_model(tmp_path, binding):
    native = Native(tmp_path)
    with pytest.raises(ValueError, match='exact requested model'):
        stage.preflight([reads(tmp_path / 'reads.fq', [FAST])], tmp_path / 'receipt.json',
                        requested=MAPPING[SUP][0], review_acknowledged=True, review_reason='Research decision',
                        reviewed_model=binding, model_mapping=MAPPING, runner=native)
    assert not native.calls
