#!/usr/bin/env python3
"""Exercise actual Medaka 2.2.0 model resolution using bundled weights; no inference."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/lib'))
import medaka_stage
import medaka
from medaka.options import basecaller_models


def main():
    default = medaka_stage.DOCUMENTED_DEFAULT['model']
    bundled = Path(medaka.__file__).parent / 'data' / (default + '_model_pt.tar.gz')
    if not bundled.is_file():
        raise RuntimeError('Expected bundled default weights are absent; this check does not provision models')
    matched = next(model for model, pair in basecaller_models.items() if pair[0] == default)
    evidence = []
    with TemporaryDirectory(prefix='hybas-medaka-selection-') as folder:
        root = Path(folder)
        for name, metadata, requested, expected in (
                ('absent', '', '', 'fallback'),
                ('unsupported', 'basecall_model_version_id=dna_r10.4.1_e8.2_400bps_fast@v4.2.0', '', 'fallback'),
                ('matched', 'basecall_model_version_id=' + matched, '', 'automatic_from_complete_read_metadata'),
                ('explicit', '', default, 'explicit')):
            reads = root / (name + '.fq')
            reads.write_text('@software_fixture ' + metadata + '\nACGT\n+\nIIII\n')
            record = medaka_stage.preflight([reads], root / (name + '.json'), requested=requested,
                selection_policy='documented_default', review_acknowledged=bool(requested),
                review_reason='Software fixture only' if requested else '', reviewed_model=requested)
            assert record['model']['selection'] == expected, record
            assert Path(record['model']['effective_path']).resolve() == bundled.resolve(), record
            assert record['model']['sha256'] == medaka_stage.sha256_file(bundled), record
            evidence.append({'case': name, 'model': record['model'],
                             'runtime_version': record['medaka_version']})
        calls = []
        def recording_runner(command, **kwargs):
            calls.append(command)
            return subprocess.run(command, **kwargs)
        try:
            medaka_stage.preflight([root / 'absent.fq'], root / 'invalid.json',
                requested='hybas_invalid_software_fixture', selection_policy='documented_default',
                review_acknowledged=True, review_reason='Deliberate invalid-model contract test',
                reviewed_model='hybas_invalid_software_fixture', runner=recording_runner)
        except subprocess.CalledProcessError:
            assert len(calls) == 1 and calls[0][-1] == 'hybas_invalid_software_fixture', calls
            assert not (root / 'invalid.json').exists()
        else:
            raise AssertionError('An invalid explicit model unexpectedly resolved')
    print(json.dumps({'checks': evidence, 'invalid_explicit_retries': 0,
        'scope': 'Actual model resolution against packaged weights only; no new model weights downloaded, '
                 'no inference, biological analysis or scientific validation.'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
