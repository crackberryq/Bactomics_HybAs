"""Dependency isolation and version attribution; no biological validation."""
import os
from pathlib import Path
import sys
import textwrap
from types import SimpleNamespace

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/lib'))
import rendering_runtime


def test_rendering_recipes_do_not_mix_incompatible_kaleido_contracts():
    qc = yaml.safe_load((ROOT / 'workflow/envs/qc.yml').read_text())['dependencies']
    report = yaml.safe_load((ROOT / 'workflow/envs/multiqc.yml').read_text())['dependencies']
    assert 'nanoplot=1.46.2' in qc and 'python-kaleido=1.0.0' in qc and 'plotly=6.1.2' in qc
    assert not any(item.startswith('multiqc') for item in qc)
    # Conda owns only the interpreter/installer. A deprecated Conda Polars
    # alias must not be mixed into the canonical pip-owned reporting stack.
    assert [item for item in report if isinstance(item, str)] == ['python=3.11', 'pip']
    pip_packages = next(item['pip'] for item in report if isinstance(item, dict))
    assert set(pip_packages) == {'multiqc==1.32', 'polars-lts-cpu==1.33.1',
                                 'plotly==6.0.1', 'kaleido==0.2.1'}


@pytest.mark.parametrize('enabled', [True, False])
def test_tool_versions_attribute_multiqc_to_its_own_manifest(tmp_path, enabled):
    content = (ROOT / 'Hybas.smk').read_text().split('rule tool_versions:', 1)[1].split('\nrule checksums:', 1)[0]
    block = textwrap.dedent(content.split('    run:\n', 1)[1].split('\n########', 1)[0])
    # A stale MultiQC row in the QC manifest must never stand in for the actual
    # reporter environment, nor imply execution when the report was disabled.
    qc, report, output = [tmp_path / name for name in ('qc.txt', 'multiqc.txt', 'versions.tsv')]
    qc.write_text('nanoplot 1.46.2 fixture\nmultiqc 0.0 stale-fixture\n')
    report.write_text('multiqc 1.32 pypi_0 pypi\n')
    inputs = {name: [] for name in ('assembly_list', 'library_list', 'pacbio_list',
                                    'preprocess_list', 'polishing_list', 'medaka_list', 'annotation_list')}
    inputs.update(qc_list=str(qc), multiqc_list=str(report) if enabled else [])
    scope = dict(os=os, input=SimpleNamespace(**inputs), output=[str(output)],
                 wildcards=SimpleNamespace(isolate='sample'), metadata_dir=lambda _: str(tmp_path),
                 base_stage=lambda _: 'unicycler')
    exec(compile(block, str(ROOT / 'Hybas.smk'), 'exec'), scope)
    rows = dict(line.split('\t') for line in output.read_text().splitlines()[1:])
    assert rows['nanoplot'] == '1.46.2'
    assert rows['multiqc'] == ('1.32' if enabled else 'NA')


def test_rendering_receipt_records_browser_identity_without_claiming_execution(tmp_path, monkeypatch):
    browser = tmp_path / 'browser'
    browser.write_bytes(b'fixture executable, never executed')
    monkeypatch.setattr(rendering_runtime, 'discover_browser', lambda: str(browser))
    monkeypatch.setattr(rendering_runtime.metadata, 'version', lambda name: 'fixture-version')
    monkeypatch.setattr(rendering_runtime.subprocess, 'run', lambda *a, **k:
                        SimpleNamespace(stdout='Fixture browser 1\n', stderr='', returncode=0))
    receipt = rendering_runtime.capture('nanoplot')
    assert receipt['browser']['resolved_path'] == str(browser)
    assert len(receipt['browser']['executable_sha256']) == 64
    assert receipt['browser']['version_output'] == 'Fixture browser 1'
    assert 'not evidence that any plot was produced' in receipt['scope']


def test_missing_browser_is_recorded_as_unavailable_without_downloading(monkeypatch):
    monkeypatch.setattr(rendering_runtime, 'discover_browser', lambda: None)
    monkeypatch.setattr(rendering_runtime.metadata, 'version', lambda name: 'fixture-version')
    def unexpected(*a, **k):
        pytest.fail('Browser was missing; no process or download should be attempted')
    monkeypatch.setattr(rendering_runtime.subprocess, 'run', unexpected)
    receipt = rendering_runtime.capture('nanoplot')
    assert receipt['browser']['status'] == 'unavailable'
