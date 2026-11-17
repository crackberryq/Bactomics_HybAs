#!/usr/bin/env bash
set -euo pipefail
# Run inside Snakemake's environment; never activate another environment here.
for tool in python medaka_consensus medaka medaka_version_report samtools minimap2 bgzip tabix bcftools; do
  location=$(command -v "$tool") || { echo "Missing executable: $tool" >&2; exit 1; }
  case "$location" in "$CONDA_PREFIX"/*) ;; *) echo "Executable outside managed environment: $tool" >&2; exit 1;; esac
done
python -m pip check
medaka --version
medaka_version_report
python - <<'PY'
from importlib import import_module, metadata
import hashlib
import json
import os
from pathlib import Path
import random
import re
import sys
from tempfile import TemporaryDirectory


def normalized(name):
    return re.sub(r'[-_.]+', '-', name).lower()


def require_inside(path, prefix):
    resolved = Path(path).resolve()
    if not resolved.is_relative_to(prefix):
        raise RuntimeError('Package file is outside the managed environment: ' + str(resolved))
    return resolved


def validate_conda_records(records):
    # Do not hide the upstream requirements mismatch by installing aliases into
    # the same Python namespaces as the canonical ONT wheels.
    aliases = {'medaka', 'mappy', 'parasail', 'parasail-python', 'ont-mappy', 'ont-parasail', 'pysam'}
    collisions = sorted({normalized(row['name']) for row in records} & aliases)
    if collisions:
        raise RuntimeError('Conflicting Conda/PyPI namespace owners: ' + ', '.join(collisions))


def verify_conda_file(path, record, prefix):
    path = require_inside(path, prefix)
    relative = path.relative_to(prefix).as_posix()
    matches = [entry for entry in record.get('paths_data', {}).get('paths', [])
               if entry.get('_path') == relative]
    if len(matches) != 1:
        raise RuntimeError('Missing unique Conda file ownership: ' + relative)
    entry = matches[0]
    # sha256_in_prefix accounts for Conda's prefix rewriting. The archive file
    # hash is usable only when no prefix replacement applies to this file.
    expected = entry.get('sha256_in_prefix') or (
        entry.get('sha256') if not entry.get('prefix_placeholder') else None)
    if not expected:
        raise RuntimeError('No installed-file checksum for Conda-owned file: ' + relative)
    with path.open('rb') as handle:
        actual = hashlib.file_digest(handle, 'sha256').hexdigest()
    if actual != expected:
        raise RuntimeError('Conda-owned file was replaced or changed: ' + relative)


def require_namespace(module, distribution, prefix):
    origin = require_inside(module.__file__, prefix)
    owners = []
    for candidate in metadata.distributions():
        if any(Path(candidate.locate_file(item)).resolve() == origin for item in candidate.files or []):
            owners.append(normalized(candidate.metadata['Name']))
    if owners != [normalized(distribution)]:
        raise RuntimeError(f'Non-unique namespace ownership for {module.__name__}: {owners}')


def main():
    prefix = Path(os.environ['CONDA_PREFIX']).resolve()
    if Path(sys.prefix).resolve() != prefix:
        raise RuntimeError('Python is outside the managed environment')
    expected = {'medaka': '2.2.0', 'ont-mappy': '2.28', 'ont-parasail': '1.3.4',
                'pysam': '0.23.0', 'torch': '2.9.1', 'numpy': '2.3.5',
                'packaging': '25.0', 'tensordict': '0.10.0'}
    actual = {name: metadata.version(name) for name in expected}
    if actual != expected:
        raise RuntimeError('Medaka package contract differs: ' + json.dumps(actual, sort_keys=True))
    aliases = {'mappy', 'parasail', 'parasail-python', 'medaka-cpu'}
    installed = {normalized(dist.metadata['Name']) for dist in metadata.distributions()}
    if aliases & installed:
        raise RuntimeError('Conflicting Python distribution aliases: ' + ', '.join(sorted(aliases & installed)))
    records = [json.loads(path.read_text()) for path in (prefix / 'conda-meta').glob('*.json')]
    validate_conda_records(records)
    modules = {name: import_module(name) for name in ('medaka', 'mappy', 'parasail', 'pysam', 'libmedaka', 'torch', 'numpy')}
    for module, distribution in (('medaka', 'medaka'), ('mappy', 'ont-mappy'), ('parasail', 'ont-parasail'),
                                 ('pysam', 'pysam'), ('libmedaka', 'medaka'), ('torch', 'torch'), ('numpy', 'numpy')):
        require_namespace(modules[module], distribution, prefix)
    # Check the actual imported entrypoints/native extensions against Conda's
    # installed payload, rather than trusting a stale package record alone.
    for package, names in (('pytorch', ('torch', 'torch._C', 'torch.version')),
                           ('numpy', ('numpy', 'numpy._core._multiarray_umath'))):
        owners = [row for row in records if row['name'] == package]
        if len(owners) != 1:
            raise RuntimeError('Missing unique Conda package owner: ' + package)
        for name in names:
            verify_conda_file(import_module(name).__file__, owners[0], prefix)

    rng = random.Random(42)
    reference = ''.join(rng.choices('ACGT', k=1600))
    aligner = modules['mappy'].Aligner(seq=reference, preset='map-ont')
    if not aligner or not any(hit.mlen >= 650 for hit in aligner.map(reference[400:1100])):
        raise RuntimeError('ONT mappy native mapping check failed')
    parasail = modules['parasail']
    if parasail.nw_trace_striped_16('ACGT', 'ACGT', 5, 2, parasail.matrix_create('ACGT', 2, -1)).score != 8:
        raise RuntimeError('ONT parasail native alignment check failed')
    with TemporaryDirectory(prefix='hybas-medaka-check-') as folder:
        fastq = Path(folder) / 'synthetic.fastq'
        fastq.write_text('@software_fixture\nACGT\n+\nIIII\n')
        with modules['pysam'].FastxFile(str(fastq)) as handle:
            read = next(handle)
            if read.sequence != 'ACGT' or read.quality != 'IIII':
                raise RuntimeError('pysam native FASTQ check failed')
    libmedaka = modules['libmedaka']
    if not hasattr(libmedaka, 'ffi') or not hasattr(libmedaka, 'lib'):
        raise RuntimeError('Medaka CFFI extension did not load')
    torch = modules['torch']
    if (torch.tensor([1.0, 2.0], device='cpu') * 2).tolist() != [2.0, 4.0]:
        raise RuntimeError('PyTorch CPU arithmetic check failed')
    print('Medaka dependency/native startup check: ' + json.dumps({
        'packages': actual, 'torch_version': torch.__version__, 'torch_cuda_build': torch.version.cuda,
        'cuda_available': torch.cuda.is_available(),
        'scope': 'Synthetic package/native-library startup only; no model inference or scientific validation. '
                 'CPU arithmetic here does not set the run device policy.',
    }, sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
PY
