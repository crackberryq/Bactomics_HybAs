"""Native CLI contracts on deterministic random reads; no biological benchmark."""
import argparse
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/lib'))
from library_inputs import preprocess_library, fastq, promote


def write_reads(path, rows):
    with gzip.open(path, 'wt') as handle:
        for name, seq in rows:
            handle.write(f'@{name}\n{seq}\n+\n' + 'I' * len(seq) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tool', choices=['flye', 'nxtrim'])
    args = parser.parse_args()
    rng = random.Random(7861)
    with tempfile.TemporaryDirectory(prefix='hybas-platform-smoke-') as tmp:
        root = Path(tmp)
        if args.tool == 'flye':
            seq = ''.join(rng.choices('ACGT', k=200000))
            path = root / 'reads.fq.gz'
            write_reads(path, [(f'read{i}', seq[(start := rng.randrange(180000)):start+20000]) for i in range(500)])
            for flag in ('--pacbio-hifi', '--pacbio-raw'):
                output = root / flag[2:]
                subprocess.run(['flye', flag, str(path), '--out-dir', str(output), '--threads', '2', '--iterations', '1'], check=True)
                promote(output / 'assembly.fasta', root / (flag[2:] + '.fasta'), 1000)
        else:
            paths = [root / 'R1.fq.gz', root / 'R2.fq.gz']
            for mate, path in enumerate(paths, 1):
                write_reads(path, [(f'pair{i}/{mate}', ''.join(rng.choices('ACGT', k=150))) for i in range(100)])
            entries, eligible, details, report, _ = preprocess_library(
                dict(library_id='nextera', library_type='mate_pair', orientation='rf', preprocessing='nextera_mate_pair'), paths, root, 2)
            assert not eligible and report.is_file()
            assert sum(details['nxtrim_read_counts'].values()) > 0
            assert all(list(fastq(p)) for entry in entries for key in ('interlaced reads', 'single reads') for p in entry.get(key, []))
        print(args.tool + ' native input/output contracts passed')


if __name__ == '__main__':
    main()
