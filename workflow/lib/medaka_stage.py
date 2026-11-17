"""Run each scheduled Medaka attempt in fresh storage with explicit model provenance.

Snakemake decides whether a stage needs running. Once scheduled, native Medaka
caches from an earlier attempt are never reused. Attempts, including failures,
remain available for diagnosis. Each file is atomically replaced; the downstream
execution recorder verifies that the receipt matches the published consensus.
"""
from __future__ import annotations
import argparse
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from library_inputs import fastq
from workflow_support import sha256_file, require_fasta, write_json

MODEL_REVIEW_CONTRACT = 2
SELECTION_POLICIES = ('strict', 'documented_default')
DOCUMENTED_DEFAULT = {
    'medaka_version': '2.2.0',
    'model': 'r1041_e82_400bps_sup_v5.2.0',
    'upstream_source_commit': '019a2ab597218a1c239fb7a5315f271a4577ed6b',
    'options_sha256': '16af07e1530bd30950cd22954e084b4d8f92316d0bc41ecd5669e768eb707f8e',
    'basis': 'upstream_non_bacterial_consensus_default',
    'metadata_match': False,
    'scientific_suitability': 'not_assessed',
}


def installed_defaults():
    from medaka import options
    return {**dict(options.default_models), 'options_sha256': sha256_file(options.__file__)}


def installed_mapping():
    try:
        from medaka.options import basecaller_models
        return dict(basecaller_models)
    except ImportError:
        return {}


def runtime_identity():
    try:
        medaka_version = version('medaka')
    except PackageNotFoundError:
        medaka_version = None
    return {'medaka_version': medaka_version,
            'device_policy': {'cuda_visible_devices': os.environ.get('CUDA_VISIBLE_DEVICES'),
                              'selection': 'native_medaka; actual inference device not observed'}}


def model_metadata(reads, known_models=()):
    models, missing, count = set(), 0, 0
    for header, _, _ in fastq(reads):
        count += 1
        found = set(re.findall(r'(?:^|\s)(?:basecall_model_version_id|basecaller_model_version_id|model_version_id)=([^\s]+)', header))
        group = re.search(r'(?:^|\s)RG:Z:([^\s]+)', header)
        if group:
            # Preserve unknown versioned IDs. Substring matching can otherwise
            # mistake an unknown version for an older unversioned model.
            explicit = re.search(r'(?:^|_)(dna_[^\s]+)', group[1])
            if explicit:
                matches = {re.sub(r'_(?:barcode[0-9]+|unclassified)$', '', explicit[1])}
            else:
                matches = {m for m in known_models if re.search(
                    r'(?:^|_)' + re.escape(m) + r'(?=$|_(?:barcode[0-9]+|unclassified)$)', group[1])}
            found.update(matches)
        if not found:
            missing += 1
        models.update(found)
    if not count:
        raise ValueError('No reads available for Medaka')
    status = 'mixed' if len(models) > 1 else 'incomplete' if missing else 'single_model_recorded'
    return {'status': status, 'inspection': 'all_reads', 'basecaller_models': sorted(models), 'reads': count,
            'reads_without_resolved_model_metadata': missing}


def model_assessment(requested, metadata, mapping):
    """Describe the installed software mapping, never scientific suitability."""
    expected = None
    if metadata['status'] != 'single_model_recorded':
        status = 'metadata_review_required'
    elif (pair := mapping.get(metadata['basecaller_models'][0])) is None or not pair[0]:
        status = 'unsupported_basecaller' if mapping else 'installed_mapping_unavailable'
    else:
        expected = pair[0]
        if not requested or requested in {expected, metadata['basecaller_models'][0] + ':consensus'}:
            status = 'matched_to_installed_mapping'
        elif requested in {model for models in mapping.values() for model in models if model}:
            status = 'mismatched_to_installed_mapping'
        else:
            status = 'unverified_explicit_model'
    return {'status': status, 'expected_consensus_model': expected,
            'basis': 'installed_medaka_basecaller_mapping', 'scientific_suitability': 'not_assessed',
            'acknowledgement_required': bool(requested and status != 'matched_to_installed_mapping')}


def review_model(requested, metadata, mapping, acknowledged, reason, reviewed_model):
    assessment = model_assessment(requested, metadata, mapping)
    assessment['review_acknowledged'] = acknowledged is True
    assessment['review_reason'] = reason
    assessment['reviewed_model'] = reviewed_model
    return assessment


def require_review(requested, metadata, assessment):
    if not requested and (metadata['status'] != 'single_model_recorded' or len(metadata['basecaller_models']) != 1):
        raise ValueError('Automatic Medaka selection requires one model identity on every read; '
                         'metadata is absent, incomplete or mixed. Review the input and set medaka_model explicitly.')
    if not requested and assessment['status'] != 'matched_to_installed_mapping':
        raise ValueError('Automatic Medaka selection is unavailable: ' + assessment['status'] +
                         '. Review an explicit model decision; no default model was substituted.')
    if assessment['acknowledgement_required']:
        if (not assessment['review_acknowledged'] or not assessment['review_reason'].strip()
                or assessment['reviewed_model'] != requested):
            raise ValueError('Explicit Medaka model requires review: ' + assessment['status'] +
                             '. Set medaka_model_review_acknowledged=true and a nonempty '
                             'medaka_model_review_reason and set medaka_model_reviewed_model to the exact requested model '
                             'in a new draft. No default model was substituted.')
        print('WARNING: Medaka explicit model review: ' + assessment['status'] +
              '; scientific suitability is not established. See model receipt.', flush=True)


def selection_decision(requested, metadata, assessment, policy, runtime):
    """Choose before native resolution; never retry a failed resolver with another model."""
    if policy not in SELECTION_POLICIES:
        raise ValueError('medaka_model_selection_policy must be strict or documented_default')
    if requested:
        require_review(requested, metadata, assessment)
        return {'selection': 'explicit', 'selected_model': requested,
                'selection_reason': 'researcher_explicit_choice'}
    if assessment['status'] == 'matched_to_installed_mapping':
        return {'selection': 'automatic_from_complete_read_metadata',
                'selected_model': assessment['expected_consensus_model'],
                'selection_reason': 'exact_installed_basecaller_mapping'}
    if policy == 'strict':
        require_review(requested, metadata, assessment)
    # Mixed identities or partial annotations cannot be described as a single
    # absent/unsupported basecaller. Require an explicit reviewed choice.
    if not metadata['basecaller_models']:
        reason = 'basecaller_metadata_absent'
    elif metadata['status'] == 'single_model_recorded' and assessment['status'] == 'unsupported_basecaller':
        reason = 'basecaller_model_unsupported'
    else:
        raise ValueError('Documented Medaka fallback is unavailable for mixed, incomplete or unavailable '
                         'model mapping evidence. Review an explicit model choice; no model was substituted.')
    if runtime.get('medaka_version') != DOCUMENTED_DEFAULT['medaka_version']:
        raise ValueError('Documented fallback requires the reviewed Medaka 2.2.0 runtime')
    # A missing/import-broken model table is a runtime failure, not absent read metadata.
    defaults = installed_defaults()
    if (defaults.get('consensus') != DOCUMENTED_DEFAULT['model']
            or defaults.get('options_sha256') != DOCUMENTED_DEFAULT['options_sha256']):
        raise ValueError('Installed Medaka default differs from the reviewed source; no fallback was substituted')
    print('WARNING: Using documented Medaka 2.2.0 consensus default ' + DOCUMENTED_DEFAULT['model'] +
          ' because ' + reason + '. This is not a metadata match or scientific validation.', flush=True)
    return {'selection': 'fallback', 'selected_model': DOCUMENTED_DEFAULT['model'],
            'selection_reason': reason, 'fallback': dict(DOCUMENTED_DEFAULT)}


def resolve_effective_model(reads, requested, metadata, runner, record=None, record_path=None, automatic_model=None,
                            decision=None):
    args = ['medaka', 'tools', 'resolve_model']
    # Resolve the exact all-read mapping decision; do not delegate to Medaka's
    # first-100-record/substr parser after this stricter inspection.
    if not requested and not automatic_model:
        raise ValueError('No checked automatic Medaka model is available')
    args += ['--model', requested or automatic_model]
    if record is not None:
        record['resolution_command'] = args
        if record_path is not None:
            write_json(record_path, record)
    result = runner(args, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    resolved = result.stdout.strip()
    if not resolved or '\n' in resolved or not Path(resolved).is_file():
        raise ValueError('Medaka did not resolve a model file; no default fallback is allowed')
    path = str(Path(resolved).resolve())
    return {'requested': requested or None, 'effective_path': path,
            'effective_model': Path(path).name, 'sha256': sha256_file(path),
            'selection': 'explicit' if requested else 'automatic_from_complete_read_metadata',
            'input_metadata': metadata, **(decision or {})}


def _selection_record(record, paths, requested, mapping, known_models, acknowledged, reason, reviewed_model):
    records = [{'path': str(path.resolve()), **model_metadata(path, known_models)} for path in paths]
    models = sorted({model for item in records for model in item['basecaller_models']})
    missing = sum(item['reads_without_resolved_model_metadata'] for item in records)
    metadata = {'status': 'mixed' if len(models) > 1 else 'incomplete' if missing else 'single_model_recorded',
                'inspection': 'all_reads', 'files': records, 'basecaller_models': models,
                'reads': sum(item['reads'] for item in records), 'reads_without_resolved_model_metadata': missing}
    record['input_metadata'] = metadata
    record['model_review'] = review_model(requested, metadata, mapping, acknowledged, reason, reviewed_model)
    # Unavailable installed mapping is an environment fault, even if read
    # annotations are also missing. Do not turn it into a default decision.
    if not mapping and not requested and record['model_selection_policy'] == 'documented_default':
        raise ValueError('Automatic Medaka selection is unavailable: installed_mapping_unavailable')
    record['model_decision'] = selection_decision(requested, metadata, record['model_review'],
                                                  record['model_selection_policy'], record)
    return metadata


def _failure(record, exc):
    record.update(status='failed', error=str(exc))
    if isinstance(exc, subprocess.CalledProcessError):
        record['failed_command'] = exc.cmd
        record['returncode'] = exc.returncode
        stderr = exc.stderr.decode(errors='replace') if isinstance(exc.stderr, bytes) else exc.stderr
        record['resolution_stderr_tail'] = (stderr or '')[-8192:]


def preflight(reads, receipt, requested='', review_acknowledged=False, review_reason='',
              runner=subprocess.run, model_mapping=None, known_models=None, reviewed_model='',
              selection_policy='strict'):
    """Inspect every raw read before costly work; execution rechecks filtered reads."""
    paths, receipt = [Path(path) for path in reads], Path(receipt)
    if not paths:
        raise ValueError('No ONT files available for Medaka preflight')
    attempts = receipt.parent / 'preflight_attempts'; attempts.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=attempts))
    record = {'schema': 'hybas-medaka-preflight-v1', 'status': 'started',
              'model_review_contract': MODEL_REVIEW_CONTRACT, 'model_selection_policy': selection_policy,
              'requested_model': requested or None, 'attempt': str(attempt.resolve()), **runtime_identity()}
    try:
        mapping = installed_mapping() if model_mapping is None else model_mapping
        known_models = tuple(mapping) if known_models is None else known_models
        record['inputs_sha256'] = {str(path.resolve()): sha256_file(path) for path in paths}
        metadata = _selection_record(record, paths, requested, mapping, known_models, review_acknowledged, review_reason, reviewed_model)
        write_json(attempt / 'execution.json', record)
        record['model'] = resolve_effective_model(paths[0], requested, metadata, runner, record, attempt / 'execution.json',
                                                  record['model_decision']['selected_model'], record['model_decision'])
        if record['inputs_sha256'] != {str(path.resolve()): sha256_file(path) for path in paths}:
            raise ValueError('ONT input changed during Medaka preflight')
        if record['model']['sha256'] != sha256_file(record['model']['effective_path']):
            raise ValueError('Medaka model changed during preflight')
        record['status'] = 'completed'
        write_json(attempt / 'execution.json', record)
        write_json(receipt, record)
        return record
    except Exception as exc:
        _failure(record, exc); write_json(attempt / 'execution.json', record)
        raise


def execute(reads, draft, output, receipt, threads, requested='', runner=subprocess.run, known_models=None,
            model_mapping=None, review_acknowledged=False, review_reason='', preflight_receipt=None, reviewed_model='',
            selection_policy='strict'):
    reads, draft, output, receipt = map(Path, (reads, draft, output, receipt))
    output.parent.mkdir(parents=True, exist_ok=True)
    attempts = output.parent / 'attempts'; attempts.mkdir(exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=attempts))
    record = {'schema': 'hybas-medaka-execution-v1', 'status': 'started',
              'model_review_contract': MODEL_REVIEW_CONTRACT, 'model_selection_policy': selection_policy,
              'attempt': str(attempt.resolve()), 'cache_policy': 'fresh_directory_per_scheduled_execution',
              'requested_model': requested or None, **runtime_identity()}
    try:
        mapping = installed_mapping() if model_mapping is None else model_mapping
        known_models = tuple(mapping) if known_models is None else known_models
        record['inputs_sha256'] = {'reads': sha256_file(reads), 'draft': sha256_file(draft)}
        metadata = _selection_record(record, [reads], requested, mapping, known_models, review_acknowledged, review_reason, reviewed_model)
        write_json(attempt / 'execution.json', record)
        record['model'] = resolve_effective_model(reads, requested, metadata, runner, record, attempt / 'execution.json',
                                                  record['model_decision']['selected_model'], record['model_decision'])
        if preflight_receipt:
            prior = json.loads(Path(preflight_receipt).read_text())
            record['preflight_receipt_sha256'] = sha256_file(preflight_receipt)
            if (prior.get('schema') != 'hybas-medaka-preflight-v1' or prior.get('status') != 'completed'
                    or prior.get('model_review_contract') != MODEL_REVIEW_CONTRACT
                    or prior.get('model_selection_policy') != selection_policy
                    or prior.get('model_decision') != record['model_decision']
                    or prior.get('requested_model') != (requested or None)
                    or prior.get('model', {}).get('sha256') != record['model']['sha256']
                    or prior.get('medaka_version') != record['medaka_version']
                    or any(prior.get('model_review', {}).get(key) != record['model_review'].get(key)
                           for key in ('review_acknowledged', 'review_reason', 'reviewed_model'))):
                raise ValueError('Medaka model or software identity differs from its preflight receipt; review a new draft')
        write_json(attempt / 'execution.json', record)
        model = record['model']
        command = ['medaka_consensus', '-t', str(threads), '-i', str(reads.resolve()),
                   '-d', str(draft.resolve()), '-o', str(attempt / 'native'), '-m', model['effective_path']]
        record['command'] = command
        write_json(attempt / 'execution.json', record)
        runner(command, check=True)
        consensus = attempt / 'native' / 'consensus.fasta'
        record['consensus'] = require_fasta(consensus)
        if record['inputs_sha256'] != {'reads': sha256_file(reads), 'draft': sha256_file(draft)}:
            raise ValueError('Medaka input changed during execution; consensus not published')
        if model['sha256'] != sha256_file(model['effective_path']):
            raise ValueError('Medaka model changed during execution; consensus not published')
        record['status'] = 'completed'
        record['consensus_sha256'] = sha256_file(consensus)
        write_json(attempt / 'execution.json', record)
        # Promotion scratch belongs to this attempt. A failed promotion must not
        # be overwritten by the next scheduled execution's scratch file.
        temporary = attempt / 'consensus.publish.tmp'
        shutil.copyfile(consensus, temporary); os.replace(temporary, output)
        receipt.parent.mkdir(parents=True, exist_ok=True)
        temporary = receipt.with_name(receipt.name + '.' + attempt.name + '.tmp')
        write_json(temporary, record); os.replace(temporary, receipt)
        return record
    except Exception as exc:
        _failure(record, exc)
        write_json(attempt / 'execution.json', record)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--reads', required=True, nargs='+')
    for key in ('draft', 'output'):
        parser.add_argument('--' + key)
    for key in ('receipt',):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--requested', default='')
    parser.add_argument('--selection-policy', choices=SELECTION_POLICIES, default='strict')
    parser.add_argument('--review-acknowledged', action='store_true')
    parser.add_argument('--review-reason', default='')
    parser.add_argument('--reviewed-model', default='')
    parser.add_argument('--preflight-receipt')
    args = vars(parser.parse_args())
    if args.pop('preflight'):
        for key in ('draft', 'output', 'threads', 'preflight_receipt'): args.pop(key)
        preflight(**args)
    else:
        if len(args['reads']) != 1 or not args['draft'] or not args['output']:
            parser.error('execution requires one reads file, --draft and --output')
        args['reads'] = args['reads'][0]
        execute(**args)


if __name__ == '__main__':
    main()
