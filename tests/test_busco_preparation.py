"""Run the sample-free database target with BUSCO's native download layout."""
import pytest
import yaml

from test_workflow import prepare, run


@pytest.fixture
def admin(tmp_path):
    config, env = prepare(tmp_path, {})
    settings = yaml.safe_load(config.read_text())
    settings.update(admin_mode=True, sample_sheet=None, busco_download_lineages=["bacteria_odb10"])
    config.write_text(yaml.safe_dump(settings))
    recorder = tmp_path / 'provider-calls'
    env['PROVIDER_CALLS'] = str(recorder)
    provider = tmp_path / "bin" / "busco"
    provider.unlink()  # prepare() installs the analysis-only stub in this fixture.
    provider.write_text("""#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

flag, lineage, option, root = sys.argv[1:]
assert (flag, lineage, option) == ('--download', 'bacteria_odb10', '--download_path')
with open(os.environ['PROVIDER_CALLS'], 'a') as log:
    print(json.dumps(sys.argv[1:]), file=log)
layout = os.environ['PROVIDER_LAYOUT']
dataset = Path(root) / ('lineages' if layout != 'flat' else '') / lineage
(dataset / 'hmms').mkdir(parents=True)
for name in ['dataset.cfg', 'hmms/example.hmm', 'lengths_cutoff', 'scores_cutoff']:
    if layout != 'incomplete' or name != 'scores_cutoff':
        (dataset / name).write_text(name)
""")
    provider.chmod(0o755)
    return env, tmp_path / 'db' / 'busco' / 'bacteria_odb10', recorder


@pytest.mark.parametrize('layout', ['lineages', 'flat'])
def test_prepare_busco_download_and_reuse(tmp_path, admin, layout):
    env, dataset, recorder = admin
    env['PROVIDER_LAYOUT'] = layout
    run(tmp_path, env, 'download_busco_db')
    for name in ['dataset.cfg', 'hmms/example.hmm', 'scores_cutoff', 'lengths_cutoff', '.complete']:
        assert (dataset / name).is_file()
    assert not list(dataset.parent.glob('.download.*'))
    # An unstamped, complete lineage must also be reusable without downloading.
    (dataset / '.complete').unlink()
    run(tmp_path, env, 'download_busco_db')
    assert (dataset / '.complete').is_file()
    assert len(recorder.read_text().splitlines()) == 1


def test_prepare_busco_rejects_incomplete_and_preserves_previous_files(tmp_path, admin):
    env, dataset, _ = admin
    dataset.mkdir(parents=True)
    previous = dataset / 'existing-file'
    previous.write_text('retain')
    env['PROVIDER_LAYOUT'] = 'incomplete'
    run(tmp_path, env, 'download_busco_db', ok=False)
    assert previous.read_text() == 'retain'
    assert not (dataset / '.complete').exists()
    assert not list(dataset.parent.glob('.download.*'))
    log = tmp_path / 'provenance/databases/busco_bacteria_odb10.log'
    assert 'scores_cutoff' in log.read_text() and 'lineages/bacteria_odb10' in log.read_text()
    env['PROVIDER_LAYOUT'] = 'lineages'
    run(tmp_path, env, 'download_busco_db')
    assert (dataset / '.complete').is_file()
    backups = list(dataset.parent.glob('.previous-incomplete.*/bacteria_odb10/existing-file'))
    assert len(backups) == 1 and backups[0].read_text() == 'retain'
