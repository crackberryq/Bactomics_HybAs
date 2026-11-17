"""Fail closed on overlapping Python owners and replaced native dependencies."""
import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def checks():
    text = (ROOT / 'workflow/envs/medaka.post-deploy.sh').read_text()
    source = text.split("python - <<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
    scope = {'__name__': 'postdeploy_test'}
    exec(compile(source, '<medaka postdeploy>', 'exec'), scope)
    return scope


def test_recipe_keeps_upstream_names_cap_and_conda_torch_constraint():
    dependencies = yaml.safe_load((ROOT / 'workflow/envs/medaka.yml').read_text())['dependencies']
    conda = [item for item in dependencies if isinstance(item, str)]
    pip = next(item['pip'] for item in dependencies if isinstance(item, dict))
    assert not {'medaka', 'mappy', 'parasail-python', 'pysam'} & {item.split('=')[0] for item in conda}
    assert {'medaka==2.2.0', 'ont-mappy==2.28', 'ont-parasail==1.3.4', 'pysam==0.23.0'} <= set(pip)
    assert 'pytorch=2.9.1' in conda and 'torch==2.9.1' in pip
    assert 'numpy=2.3.5' in conda and 'numpy==2.3.5' in pip
    assert not any('cpu' in item or 'cuda' in item for item in conda + pip)


@pytest.mark.parametrize('alias', ['mappy', 'parasail-python', 'medaka', 'pysam'])
def test_conda_alias_rejected_even_if_pip_check_would_succeed(alias):
    with pytest.raises(RuntimeError, match='namespace owners'):
        checks()['validate_conda_records']([{'name': alias}])


def test_duplicate_distribution_owner_rejected(tmp_path, monkeypatch):
    scope = checks()
    origin = tmp_path / 'mappy.so'
    origin.write_bytes(b'fixture only')
    def dist(name):
        return SimpleNamespace(metadata={'Name': name}, files=[origin.name], locate_file=lambda _: origin)
    monkeypatch.setattr(scope['metadata'], 'distributions', lambda: [dist('ont-mappy'), dist('mappy')])
    with pytest.raises(RuntimeError, match='Non-unique namespace'):
        scope['require_namespace'](SimpleNamespace(__file__=origin, __name__='mappy'), 'ont-mappy', tmp_path)
    monkeypatch.setattr(scope['metadata'], 'distributions', lambda: [dist('ont-mappy')])
    scope['require_namespace'](SimpleNamespace(__file__=origin, __name__='mappy'), 'ont-mappy', tmp_path)


def test_installed_conda_checksum_detects_replacement_and_honors_prefix_rewrite(tmp_path):
    verify = checks()['verify_conda_file']
    origin = tmp_path / 'torch.so'
    origin.write_bytes(b'installed prefix bytes')
    entry = {'_path': origin.name, 'sha256': hashlib.sha256(b'archive bytes').hexdigest(),
             'prefix_placeholder': '/old/prefix',
             'sha256_in_prefix': hashlib.sha256(origin.read_bytes()).hexdigest()}
    record = {'paths_data': {'paths': [entry]}}
    verify(origin, record, tmp_path)
    origin.write_bytes(b'replaced by another installer')
    with pytest.raises(RuntimeError, match='replaced or changed'):
        verify(origin, record, tmp_path)
    del entry['sha256_in_prefix']
    with pytest.raises(RuntimeError, match='No installed-file checksum'):
        verify(origin, record, tmp_path)


def test_unowned_files_and_symlink_escape_are_rejected(tmp_path):
    verify = checks()['verify_conda_file']
    prefix = tmp_path / 'env'
    prefix.mkdir()
    inside = prefix / 'torch.so'
    inside.write_bytes(b'fixture')
    with pytest.raises(RuntimeError, match='unique Conda file ownership'):
        verify(inside, {'paths_data': {'paths': []}}, prefix)
    outside = tmp_path / 'outside.so'
    outside.write_bytes(b'fixture')
    inside.unlink()
    inside.symlink_to(outside)
    with pytest.raises(RuntimeError, match='outside the managed environment'):
        verify(inside, {'paths_data': {'paths': []}}, prefix)
