"""Missing-browser authorization tests; actual downloads/rendering run in CI."""
import ast
from importlib import metadata
from pathlib import Path
import subprocess
import sys
from types import ModuleType, SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def preparation_block():
    text = (ROOT / 'workflow/envs/qc.post-deploy.sh').read_text().split("python - <<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]
    module = ast.parse(text)
    # Execute the shipped dependency/browser preparation, ending before the
    # actual rendering fixture. No duplicated preparation implementation.
    stop = next(index for index, node in enumerate(module.body) if isinstance(node, ast.With))
    return compile(ast.Module(body=module.body[:stop], type_ignores=[]), '<native browser preparation>', 'exec')


@pytest.mark.parametrize('allowed,existing,downloads', [('0', None, 0), ('1', None, 1), ('1', '/fixture/chrome', 0)])
def test_browser_preparation_requires_explicit_opt_in_and_reuses_existing(monkeypatch, allowed, existing, downloads):
    monkeypatch.setenv('HYBAS_PREPARE_BROWSER', allowed)
    monkeypatch.delenv('BROWSER_PATH', raising=False)
    versions = {'nanoplot': '1.46.2', 'plotly': '6.1.2', 'kaleido': '1.0.0'}
    monkeypatch.setattr(metadata, 'version', versions.__getitem__)
    chromium = ModuleType('choreographer.browsers.chromium')
    chromium.Chromium = SimpleNamespace(find_browser=lambda **kwargs: existing)
    monkeypatch.setitem(sys.modules, 'choreographer.browsers.chromium', chromium)
    calls = []
    kaleido = ModuleType('kaleido')
    def provision():
        calls.append('explicit-download')
        return '/fixture/downloaded-chrome'
    kaleido.get_chrome_sync = provision
    monkeypatch.setitem(sys.modules, 'kaleido', kaleido)
    monkeypatch.setattr(subprocess, 'check_output', lambda *a, **kw: 'Fixture Chrome version\n')
    if allowed == '0' and existing is None:
        with pytest.raises(RuntimeError, match='No browser is downloaded during checks'):
            exec(preparation_block(), {})
    else:
        exec(preparation_block(), {})
    assert len(calls) == downloads
