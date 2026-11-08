"""Content identity for an immutable, prepared BUSCO lineage directory."""
import argparse
import hashlib
import json
from pathlib import Path
from workflow_support import sha256_file, write_json


def inventory(path):
    root = Path(path)
    files = []
    for child in sorted(root.rglob('*')):
        if child.is_file() and child.name != '.complete':
            files.append({'path': child.relative_to(root).as_posix(), 'sha256': sha256_file(child)})
    body = json.dumps(files, sort_keys=True, separators=(',', ':')).encode()
    return {'schema': 'hybas-busco-dataset-v1', 'directory': str(root.resolve()),
            'files': files, 'sha256': hashlib.sha256(body).hexdigest()}


def verify(path, expected, receipt=None):
    current = inventory(path)
    if not current['files'] or current['sha256'] != expected:
        raise ValueError('BUSCO lineage contents changed after workflow planning; rerun with a stable prepared dataset')
    if receipt:
        write_json(receipt, current)
    return current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--path', required=True)
    parser.add_argument('--expected', required=True)
    parser.add_argument('--receipt')
    verify(**vars(parser.parse_args()))


if __name__ == '__main__':
    main()
