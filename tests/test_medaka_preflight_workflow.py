"""The real Snakemake graph stops before expensive tools on model review failures."""
from pathlib import Path

import pytest
import yaml

from test_workflow import prepare, run, calls


@pytest.mark.parametrize('requested', ['', 'explicit-research-model'])
def test_model_review_blocks_before_read_processing_and_assembly(tmp_path, requested):
    config, env = prepare(tmp_path, {'sample': 'hybrid'}, kraken=True)
    cfg = yaml.safe_load(config.read_text())
    cfg.update(medaka_model=requested, medaka_model_review_acknowledged=False,
               medaka_model_review_reason='')
    config.write_text(yaml.safe_dump(cfg))
    run(tmp_path, env, ok=False)
    expensive = {'pigz', 'NanoPlot', 'kraken2', 'unicycler', 'spades.py', 'racon', 'medaka_consensus'}
    assert not expensive.intersection(call['tool'] for call in calls(tmp_path))
    records = list((tmp_path / 'sample/work/polishing/medaka/preflight_attempts').glob('*/execution.json'))
    assert len(records) == 1
    text = (tmp_path / 'sample/logs/medaka_model_preflight.log').read_text()
    assert ('requires review' if requested else 'requires one model identity') in text
