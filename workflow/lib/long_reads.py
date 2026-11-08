"""Platform-labelled long-read transport and counts; no taxonomy policy."""
import argparse
import gzip
import json
from pathlib import Path

from library_inputs import fastq, record, write_json


def prepare(plan, platform, output, audit):
    audit = Path(audit)
    audit.mkdir(parents=True, exist_ok=True)
    counts = []
    with gzip.open(output, 'wt') as reads, gzip.open(audit / 'original_headers.jsonl.gz', 'wt') as ids:
        for li, lib in enumerate(plan['libraries']):
            if lib['platform'] != platform:
                continue
            for gi, group in enumerate(lib['read_groups']):
                total, bases = 0, 0
                for total, value in enumerate(fastq(group['reads']), 1):
                    name = f'hybas_{platform}_L{li}G{gi}R{total}'
                    record(reads, value, name)
                    ids.write(json.dumps({'id': name, 'original_header': value[0]}) + '\n')
                    bases += len(value[1])
                if not total:
                    raise ValueError(f"{lib['library_id']}/{group['read_group_id']}: no {platform} reads")
                counts.append({'library_id': lib['library_id'], 'read_group_id': group['read_group_id'],
                               'reads': total, 'bases': bases, 'namespace': f'hybas_{platform}_L{li}G{gi}'})
    write_json(audit / 'plan.json', plan)
    write_json(audit / 'read_groups.json', counts)


def metrics(raw, cleaned, chemistry, output):
    def count(path):
        reads = bases = 0
        for reads, (_, seq, _) in enumerate(fastq(path), 1):
            bases += len(seq)
        return {'reads': reads, 'bases': bases}
    raw_counts, clean_counts = count(raw), count(cleaned)
    if not clean_counts['reads']:
        raise ValueError('No PacBio reads remain after cleaning')
    write_json(output, {'platform': 'pacbio', 'chemistry': chemistry,
                        'raw': raw_counts, 'cleaned': clean_counts,
                        'preprocessing': 'full_fastq_audit_then_AWARE; no_ONT_length_or_quality_filter'})


def dataset(source, output, ont):
    data = json.loads(Path(source).read_text())
    if ont:
        data.append({'type': 'nanopore', 'single reads': [str(Path(p).absolute()) for p in ont]})
    write_json(output, data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    for key in ('plan', 'platform', 'output', 'audit'):
        prep.add_argument('--' + key, required=True)
    count = sub.add_parser('metrics')
    for key in ('raw', 'cleaned', 'chemistry', 'output'):
        count.add_argument('--' + key, required=True)
    ds = sub.add_parser('dataset')
    ds.add_argument('--source', required=True)
    ds.add_argument('--output', required=True)
    ds.add_argument('--ont', nargs='*', default=[])
    args = vars(parser.parse_args()); action = args.pop('action')
    if action == 'prepare':
        args['plan'] = json.loads(args['plan'])
        prepare(**args)
    elif action == 'metrics':
        metrics(**args)
    else:
        dataset(**args)


if __name__ == '__main__':
    main()
