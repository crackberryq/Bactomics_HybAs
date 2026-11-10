"""File-preservation contracts using synthetic data; no biological tools run."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
import medaka_stage


def test_medaka_failed_publication_scratch_survives_later_attempt(tmp_path, monkeypatch):
    reads, draft, model = [tmp_path / name for name in ('reads.fq', 'draft.fa', 'model')]
    reads.write_text('@fixture model_version_id=model\nACGT\n+\nIIII\n')
    draft.write_text('>draft\nACGT\n')
    model.write_text('synthetic model')
    output, receipt = tmp_path / 'consensus.fa', tmp_path / 'receipt.json'
    old_scratch = tmp_path / 'consensus.fa.tmp'
    old_scratch.write_text('previous incomplete publication')

    def runner(args, **kwargs):
        if args[0] == 'medaka':
            return SimpleNamespace(stdout=str(model) + '\n')
        native = Path(args[args.index('-o') + 1])
        native.mkdir()
        (native / 'consensus.fasta').write_text('>fixture\nACGT\n')
        return SimpleNamespace(returncode=0)

    replace = medaka_stage.os.replace

    def fail_publication(source, destination):
        if Path(destination) == output:
            raise OSError('synthetic publication failure')
        return replace(source, destination)

    with monkeypatch.context() as patch:
        patch.setattr(medaka_stage.os, 'replace', fail_publication)
        with pytest.raises(OSError, match='synthetic publication failure'):
            medaka_stage.execute(reads, draft, output, receipt, 1, runner=runner, known_models=())
    failed = next((tmp_path / 'attempts').iterdir())
    saved = (failed / 'consensus.publish.tmp').read_bytes()
    assert saved == (failed / 'native/consensus.fasta').read_bytes()
    assert json.loads((failed / 'execution.json').read_text())['status'] == 'failed'
    result = medaka_stage.execute(reads, draft, output, receipt, 1, runner=runner, known_models=())
    assert Path(result['attempt']) != failed
    assert (failed / 'consensus.publish.tmp').read_bytes() == saved
    assert old_scratch.read_text() == 'previous incomplete publication'
    assert output.read_bytes() == saved

