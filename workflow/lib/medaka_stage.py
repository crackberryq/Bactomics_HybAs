"""Run each scheduled Medaka attempt in fresh storage with explicit model provenance.

Snakemake decides whether a stage needs running. Once scheduled, native Medaka
caches from an earlier attempt are never reused. Attempts, including failures,
remain available for diagnosis. Each file is atomically replaced; the downstream
execution recorder verifies that the receipt matches the published consensus.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from library_inputs import fastq
from workflow_support import sha256_file, require_fasta, write_json


def model_metadata(reads, known_models=()):
    models, missing, count = set(), 0, 0
    for header, _, _ in fastq(reads):
        count += 1
        found = set(re.findall(r'(?:^|\s)(?:basecall_model_version_id|model_version_id)=([^\s]+)', header))
        group = re.search(r'(?:^|\s)RG:Z:([^\s]+)', header)
        if group:
            matches = {m for m in known_models if m in group[1]}
            # A versioned model can contain an older unversioned model name.
            matches = {m for m in matches if not any(m != n and m in n for n in matches)}
            found.update(matches)
        if not found:
            missing += 1
        models.update(found)
    if not count:
        raise ValueError('No reads available for Medaka')
    status = 'mixed' if len(models) > 1 else 'incomplete' if missing else 'single_model_recorded'
    return {'status': status, 'basecaller_models': sorted(models), 'reads': count,
            'reads_without_resolved_model_metadata': missing}


def resolve_effective_model(reads, requested, metadata, runner):
    if not requested and (metadata['status'] != 'single_model_recorded' or len(metadata['basecaller_models']) != 1):
        raise ValueError('Automatic Medaka selection requires one model identity on every read; '
                         'metadata is absent, incomplete or mixed. Review the input and set medaka_model explicitly.')
    args = ['medaka', 'tools', 'resolve_model']
    args += ['--model', requested] if requested else ['--auto_model', 'consensus', str(reads)]
    result = runner(args, check=True, text=True, stdout=subprocess.PIPE)
    resolved = result.stdout.strip()
    if not resolved or '\n' in resolved or not Path(resolved).is_file():
        raise ValueError('Medaka did not resolve a model file; no default fallback is allowed')
    path = str(Path(resolved).resolve())
    return {'requested': requested or None, 'effective_path': path,
            'effective_model': Path(path).name, 'sha256': sha256_file(path),
            'selection': 'explicit' if requested else 'automatic_from_complete_read_metadata',
            'input_metadata': metadata}


def execute(reads, draft, output, receipt, threads, requested='', runner=subprocess.run, known_models=None):
    reads, draft, output, receipt = map(Path, (reads, draft, output, receipt))
    output.parent.mkdir(parents=True, exist_ok=True)
    attempts = output.parent / 'attempts'; attempts.mkdir(exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix='attempt-', dir=attempts))
    record = {'schema': 'hybas-medaka-execution-v1', 'status': 'started',
              'attempt': str(attempt.resolve()), 'cache_policy': 'fresh_directory_per_scheduled_execution'}
    try:
        if known_models is None:
            try:
                from medaka.options import basecaller_models
                known_models = tuple(basecaller_models)
            except ImportError:
                known_models = ()
        record['inputs_sha256'] = {'reads': sha256_file(reads), 'draft': sha256_file(draft)}
        metadata = model_metadata(reads, known_models)
        record['model'] = resolve_effective_model(reads, requested, metadata, runner)
        if requested and metadata['status'] != 'single_model_recorded':
            print('Medaka: explicit model selected; read metadata is ' + metadata['status'] + '. See execution receipt.', flush=True)
        write_json(attempt / 'execution.json', record)
        model = record['model']
        command = ['medaka_consensus', '-t', str(threads), '-i', str(reads.resolve()),
                   '-d', str(draft.resolve()), '-o', str(attempt / 'native'), '-m', model['effective_path']]
        record['command'] = command
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
        temporary = output.with_name(output.name + '.tmp')
        shutil.copyfile(consensus, temporary); os.replace(temporary, output)
        receipt.parent.mkdir(parents=True, exist_ok=True)
        temporary = receipt.with_name(receipt.name + '.tmp')
        write_json(temporary, record); os.replace(temporary, receipt)
        return record
    except Exception as exc:
        record.update(status='failed', error=str(exc))
        write_json(attempt / 'execution.json', record)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('reads', 'draft', 'output', 'receipt'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--threads', type=int, required=True)
    parser.add_argument('--requested', default='')
    execute(**vars(parser.parse_args()))


if __name__ == '__main__':
    main()
