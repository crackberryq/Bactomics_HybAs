"""Synthetic contracts for library identity, cleaning, assembly and resume."""
import gzip
import json
from pathlib import Path
import sys

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
from input_manifest import load_manifest
from library_inputs import execution_samples, pairs, prepare as prepare_reads, trim
from test_workflow import prepare, run, calls, computational_runs, settings_snapshots


def manifest_run(tmp, strategy='dual_illumina_complementary', polish='tripolish', kraken=False):
    cfg, env = prepare(tmp, {'sample': 'illumina_only'}, polish=polish, kraken=kraken)
    libs = []
    for li, size in enumerate((600, 250)):
        groups = []
        for gi in range(2 if li == 0 else 1):
            files = []
            for mate in (1, 2):
                # Arbitrary filenames and duplicate original IDs across all groups.
                path = tmp / 'reads with spaces' / f'lib{li}_group{gi}_mate{mate}.fq.gz'
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(gzip.compress(f'@duplicate/{mate}\nACGT\n+\nIIII\n'.encode()))
                files.append(str(path))
            groups.append(dict(read_group_id=f'group{gi}', run_id=f'run{gi}', r1=files[0], r2=files[1]))
        libs.append(dict(library_id=f'lib{li}', platform='illumina', library_type='paired_end',
                         orientation='fr', insert_size_bp=size, read_groups=groups))
    if strategy == 'illumina_only':
        libs = [libs[1]]
    if strategy == 'repeat_one_library':
        libs, strategy = [libs[0]], 'dual_illumina_runs'
    manifest = tmp / 'inputs.yaml'
    manifest.write_text(yaml.safe_dump(dict(schema_version=1, samples=[dict(sample_id='sample', mode=strategy, libraries=libs)])))
    config = yaml.safe_load(cfg.read_text())
    config.update(input_manifest=str(manifest), sample_sheet='')
    cfg.write_text(yaml.safe_dump(config))
    return cfg, env, manifest


@pytest.mark.parametrize('strategy,polish,count', [
    ('dual_illumina_complementary', 'tripolish', 2),
    ('dual_illumina_runs', 'none', 2),
    ('repeat_one_library', 'none', 1),
    ('illumina_only', 'none', 1),
])
def test_library_routes_and_noop(tmp_path, strategy, polish, count):
    _, env, _ = manifest_run(tmp_path, strategy, polish)
    run(tmp_path, env)
    assert not any(c['tool'] == 'unicycler' for c in calls(tmp_path))
    assert sum(c['tool'] == 'fastp' for c in calls(tmp_path)) == count
    assert sum(c['tool'] == 'spades.py' for c in calls(tmp_path)) == 1
    bundle = tmp_path / 'sample/trimmed/libraries'
    dataset = yaml.safe_load((bundle / 'dataset.yaml').read_text())
    metadata = json.loads((bundle / 'libraries.json').read_text())
    assert len(dataset) == count
    if count == 2:
        assert metadata['spades_library_order'] == ['lib1', 'lib0']
        assert [x['trimmed_pairs'] for x in metadata['libraries']] == [2, 1]
    audit = tmp_path / 'sample/work/inputs_merged/library_identity'
    with gzip.open(audit / 'original_headers.jsonl.gz', 'rt') as handle:
        headers = [json.loads(x) for x in handle]
    assert len({x['id'] for x in headers}) == len(headers)
    assert {x['r1'] for x in headers} == {'@duplicate/1'}
    execution = json.loads((tmp_path / 'sample/metadata/polish_execution.json').read_text())
    assert execution['stages'][0]['stage'] == 'spades'
    assert execution['final_source_stage'] == ('polypolish' if polish == 'tripolish' else 'spades')
    final = (tmp_path / 'sample/work/final/assembly.final.fasta').read_text()
    assert '>short' not in final
    assert final.strip().endswith('CT' if polish == 'tripolish' else 'CC')
    summary = json.loads((tmp_path / 'sample/metadata/run_summary.json').read_text())
    assert summary['assembler'] == 'spades' and summary['library_count'] == count
    frozen = json.loads(settings_snapshots(tmp_path)[0].read_text())
    assert frozen['sample_plans']['sample']['library_plan']['assembler'] == 'spades'
    assert frozen['request_context']['input_manifest_sha256']
    jobs = computational_runs(tmp_path)[0]['jobs']
    assert {'merge_illumina_libraries', 'fastp_libraries', 'spades_libraries'} <= {x['rule'] for x in jobs}
    before = calls(tmp_path)
    run(tmp_path, env)
    assert calls(tmp_path) == before
    assert computational_runs(tmp_path)[-1]['status'] == 'no_work'
    if strategy == 'dual_illumina_complementary':
        history = computational_runs(tmp_path)
        run(tmp_path, env, '--report', str(tmp_path / 'report.zip'))
        assert computational_runs(tmp_path) == history
        from zipfile import ZipFile
        with ZipFile(tmp_path / 'report.zip') as archive:
            html = [archive.read(x).decode() for x in archive.namelist() if x.endswith('.html')]
        assert any('Illumina library trimming' in x and 'lib0' in x and 'lib1' in x for x in html)


def test_library_failure_resume_and_polish_change_reuse(tmp_path):
    cfg, env, _ = manifest_run(tmp_path, polish='none', kraken=True)
    run(tmp_path, dict(env, HYBAS_TEST_FAIL_TOOL='spades.py'), ok=False)
    before = calls(tmp_path)
    run(tmp_path, env)
    new = calls(tmp_path)[len(before):]
    assert any(c['tool'] == 'spades.py' for c in new)
    assert not any(c['tool'] in {'fastp', 'kraken2'} for c in new)
    before = calls(tmp_path)
    config = yaml.safe_load(cfg.read_text()); config['polish_mode'] = 'tripolish'; cfg.write_text(yaml.safe_dump(config))
    run(tmp_path, env)
    new = calls(tmp_path)[len(before):]
    assert any(c['tool'] == 'polypolish' for c in new)
    assert not any(c['tool'] in {'fastp', 'kraken2', 'spades.py'} for c in new)


def test_manifest_setup_preview_selects_library_environment(tmp_path):
    _, env, _ = manifest_run(tmp_path)
    result = run(tmp_path, env, 'setup_envs', '--config', 'setup_only=true', '-n')
    assert 'capture_env_library_assembly_global' in result.stdout
    assert calls(tmp_path) == [] and not settings_snapshots(tmp_path)


@pytest.mark.parametrize('change', ['mate_pair', 'rf', 'ff', 'pacbio', 'sheet_mismatch'])
def test_unsupported_manifest_fails_before_jobs(tmp_path, change):
    cfg, env, manifest = manifest_run(tmp_path)
    data = yaml.safe_load(manifest.read_text())
    lib = data['samples'][0]['libraries'][0]
    if change == 'mate_pair':
        lib['library_type'], lib['orientation'] = 'mate_pair', 'rf'
    elif change in {'rf', 'ff'}:
        lib['orientation'] = change
    elif change == 'pacbio':
        data['samples'][0].update(mode='pacbio_only', libraries=[dict(library_id='hifi', platform='pacbio',
            chemistry='hifi', read_groups=[dict(read_group_id='group', run_id='run', reads=lib['read_groups'][0]['r1'])])])
    else:
        sheet = tmp_path / 'other.tsv'; sheet.write_text('isolate\nother\n')
        config = yaml.safe_load(cfg.read_text()); config['sample_sheet'] = str(sheet); cfg.write_text(yaml.safe_dump(config))
    manifest.write_text(yaml.safe_dump(data))
    result = run(tmp_path, env, ok=False)
    assert 'metadata-only' in result.stdout or 'exactly the same sample IDs' in result.stdout
    assert calls(tmp_path) == [] and not settings_snapshots(tmp_path)


@pytest.mark.parametrize('bad', ['mismatch', 'wrong_mate', 'truncated', 'unequal'])
def test_full_pair_audit(tmp_path, bad):
    _, _, manifest = manifest_run(tmp_path)
    plan = execution_samples(load_manifest(manifest))['sample']
    path = Path(plan['libraries'][0]['read_groups'][0]['r2'])
    prefix = '@duplicate/2\nACGT\n+\nIIII\n'
    extra = {'mismatch': '@other/2\nACGT\n+\nIIII\n', 'wrong_mate': '@duplicate/1\nACGT\n+\nIIII\n',
             'truncated': '@duplicate/2\nACGT\n+\nI\n', 'unequal': ''}[bad]
    path.write_bytes(gzip.compress((prefix + extra).encode()))
    left = Path(plan['libraries'][0]['read_groups'][0]['r1'])
    left.write_bytes(gzip.compress(b'@duplicate/1\nACGT\n+\nIIII\n' * 2))
    with pytest.raises(ValueError):
        prepare_reads(plan, tmp_path / 'r1.gz', tmp_path / 'r2.gz', tmp_path / 'audit')


def test_cleaning_restores_libraries_and_rejects_empty_library(tmp_path):
    _, env, manifest = manifest_run(tmp_path)
    plan = execution_samples(load_manifest(manifest))['sample']
    r1, r2 = tmp_path / 'merged1.gz', tmp_path / 'merged2.gz'
    prepare_reads(plan, r1, r2, tmp_path / 'audit')
    # Retain one pair from each preparation after a synthetic cleaning decision.
    selected = list(pairs(r1, r2))[1:]
    for mate, path in enumerate((r1, r2)):
        with gzip.open(path, 'wt') as out:
            for pair in selected:
                header, seq, qual = pair[mate]; out.write(f'{header}\n{seq}\n+\n{qual}\n')
    import os
    from unittest.mock import patch
    with patch.dict(os.environ, env):
        trim(plan, r1, r2, tmp_path / 'trim', tmp_path / 'out1.gz', tmp_path / 'out2.gz',
             tmp_path / 'report.html', tmp_path / 'report.json', 1)
    assert [x['cleaned_pairs'] for x in json.loads((tmp_path / 'report.json').read_text())['libraries']] == [1, 1]
    for mate, path in enumerate((r1, r2)):
        with gzip.open(path, 'wt') as out:
            header, seq, qual = selected[-1][mate]; out.write(f'{header}\n{seq}\n+\n{qual}\n')
    with pytest.raises(ValueError, match='no pairs remain after cleaning'):
        trim(plan, r1, r2, tmp_path / 'empty', tmp_path / 'out1.gz', tmp_path / 'out2.gz',
             tmp_path / 'report.html', tmp_path / 'report.json', 1)
