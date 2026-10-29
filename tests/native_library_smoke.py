"""Native SPAdes dataset/contig contract using deterministic random toy reads.

No organism, reference genome or biological accuracy claim is involved.
Run in workflow/envs/library_assembly.yml; outputs live in a temporary directory.
"""
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
from library_inputs import promote


def main():
    rng = random.Random(49013)
    sequence = ''.join(rng.choices('ACGT', k=10000))
    with tempfile.TemporaryDirectory(prefix='hybas-library-smoke-') as temp:
        root = Path(temp)
        dataset = []
        for li, insert in enumerate((300, 600)):
            r1, r2 = root / f'lib{li}_1.fq.gz', root / f'lib{li}_2.fq.gz'
            with gzip.open(r1, 'wt') as left, gzip.open(r2, 'wt') as right:
                for index in range(1000):
                    start = rng.randrange(len(sequence) - insert)
                    a = sequence[start:start + 150]
                    b = sequence[start + insert - 150:start + insert].translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                    for handle, seq in ((left, a), (right, b)):
                        handle.write(f'@lib{li}_pair{index}\n{seq}\n+\n'+ 'I'*150 + '\n')
            dataset.append({'type': 'paired-end', 'orientation': 'fr',
                            'left reads': [str(r1)], 'right reads': [str(r2)]})
        source = root / 'dataset.yaml'; source.write_text(json.dumps(dataset))
        subprocess.run(['spades.py', '--dataset', str(source), '-o', str(root / 'assembly'),
                        '-t', '2', '-m', '2'], check=True)
        promote(root / 'assembly/contigs.fasta', root / 'stage_spades.fasta', 1000)
        assert (root / 'stage_spades.fasta').stat().st_size > 1000
        print('Native multi-library dataset and contig promotion contracts passed')


if __name__ == '__main__':
    main()
