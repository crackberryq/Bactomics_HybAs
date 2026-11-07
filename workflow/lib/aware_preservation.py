"""Experimental, lossless taxonomic read preservation; no strain reconstruction.

Consumes a HybAs identification-review receipt. Classification labels are evidence
categories, never proof that reads originate from a single genome.
"""
from collections import Counter, OrderedDict
from datetime import datetime, timezone
import csv
import gzip
import html
import itertools
import json
from pathlib import Path
import platform as runtime_platform
import re
import resource
import shutil
import sqlite3
import tempfile
import time

from identification_review import SCHEMA as REVIEW_SCHEMA, digest, identity, read_json

VERSION = 'aware-preservation-0.2.0a1'
SCHEMA = 'hybas-aware-plus-preservation-v1'
LIMITATION = ('Taxonomic assignments are provisional. Species bins can contain multiple strains '
              'or misclassified reads. Long-read linkage has not been evaluated. These are read '
              'groups, not reconstructed isolates or guaranteed genomes. Ambiguous reads are '
              'stored once; candidate associations are incomplete when taxonomy is unresolved.')
FILES = {'illumina_report': 'illumina.report.txt', 'ont_report': 'nanopore.report.txt',
         'illumina_calls': 'illumina.kraken2.txt', 'ont_calls': 'nanopore.kraken2.txt',
         'r1': 'illumina_R1.merged.fq.gz', 'r2': 'illumina_R2.merged.fq.gz',
         'ont': 'ont_merged.fastq.gz'}


def taxonomy(reports):
    """Use the report's indentation, not rank letters, to retain ancestry."""
    nodes, counts = {}, {}
    for platform, path in reports.items():
        stack, direct, seen = [], Counter(), set()
        with Path(path).open() as stream:
            for line in stream:
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) != 6:
                    raise ValueError('Expected the standard six-column Kraken report')
                pct, clade, assigned, rank, taxid, label = parts
                if not taxid.isdigit() or taxid in seen or not assigned.isdigit() or not clade.isdigit():
                    raise ValueError('Invalid or repeated taxon in Kraken report')
                seen.add(taxid); direct[taxid] = int(assigned)
                if taxid == '0':
                    continue
                depth = len(label) - len(label.lstrip(' '))
                while stack and stack[-1][0] >= depth:
                    stack.pop()
                parent = stack[-1][1] if stack else None
                node = {'taxid': taxid, 'name': label.strip(), 'rank': rank, 'parent': parent}
                if not node['name'] or (parent is None and taxid != '1'):
                    raise ValueError('Kraken report is missing the rooted taxonomy hierarchy')
                if taxid in nodes and nodes[taxid] != node:
                    raise ValueError('Kraken reports disagree on taxon ancestry or identity')
                nodes[taxid] = node; stack.append((depth, taxid))
        counts[platform] = +direct
    return nodes, counts


def lineage(nodes, taxid):
    chain = []
    while taxid is not None:
        if taxid not in nodes or taxid in chain:
            raise ValueError('Incomplete or cyclic report taxonomy')
        chain.append(taxid); taxid = nodes[taxid]['parent']
    return chain


def route(nodes, species_for, descendants, taxid, hitlist):
    tokens = hitlist.split()
    if not tokens:
        raise ValueError('Empty Kraken hit list')
    for token in tokens:
        if token == '|:|':
            continue
        match = re.fullmatch(r'(\d+|A):(\d+)', token)
        if not match or int(match[2]) < 1:
            raise ValueError('Unsupported Kraken nucleotide hit list')
    if taxid == '0':
        return '_unclassified', 'unclassified', [], False
    if taxid not in nodes:
        raise ValueError('Assignment TaxID is absent from its reports')
    selected = species_for.get(taxid)
    candidates = set(descendants.get(taxid, ()))
    observed, unresolved, informative = set(), False, False
    for token in tokens:
        if token == '|:|':
            continue
        match = re.fullmatch(r'(\d+|A):(\d+)', token)
        hit = match[1]
        if hit in {'0', 'A'}:
            continue
        informative = True
        if hit not in nodes:
            unresolved = True
            continue
        if hit in species_for:
            observed.add(species_for[hit])
        elif selected is not None and hit not in lineage(nodes, selected):
            # A conflicting higher-rank branch is not evidence for the selected species.
            unresolved = True
            candidates.update(descendants.get(hit, ()))
    candidates.update(observed)
    if selected:
        candidates.add(selected)
    if unresolved:
        return '_ambiguous', 'unresolved_hit_taxonomy', sorted(candidates), True
    if selected and informative and observed <= {selected}:
        return selected, 'kraken_species_assignment', [selected], False
    reason = 'conflicting_species_hits' if selected else 'assignment_above_species'
    if not informative:
        reason = 'no_informative_hits'
    return '_ambiguous', reason, sorted(candidates), False


def read_key(value, paired=False):
    key = value.split()[0] if value.split() else ''
    if paired and key.endswith(('/1', '/2')):
        key = key[:-2]
    if not key or any(c in key for c in '\t\r\n\x00'):
        raise ValueError('Missing or invalid read identifier')
    return key


def fastq(path):
    """Read the four-line FASTQ produced by the current HybAs merge stage."""
    if Path(path).stat().st_size == 0:
        return
    with gzip.open(path, 'rb') as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            sequence, plus, quality = (stream.readline() for _ in range(3))
            if (not header.startswith(b'@') or not plus.startswith(b'+') or not quality
                    or not sequence.rstrip(b'\r\n')
                    or len(sequence.rstrip(b'\r\n')) != len(quality.rstrip(b'\r\n'))):
                raise ValueError('Malformed or truncated four-line FASTQ')
            name = header[1:].decode('utf-8').strip()
            if plus[1:].strip() and plus[1:].strip() != header[1:].strip():
                raise ValueError('FASTQ separator identifier differs from its header')
            yield name, header + sequence + plus + quality, len(sequence.rstrip(b'\r\n'))


class Writers:
    """Limit open files; concatenated gzip members preserve valid FASTQ streams."""
    def __init__(self, root, limit=32):
        self.root, self.limit, self.opened = root, limit, OrderedDict()

    def write(self, relative, data):
        if relative not in self.opened:
            if len(self.opened) >= self.limit:
                _, pair = self.opened.popitem(last=False)
                pair[0].close(); pair[1].close()
            path = self.root / relative; path.parent.mkdir(parents=True, exist_ok=True)
            raw = path.open('ab')
            self.opened[relative] = (gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0), raw)
        self.opened.move_to_end(relative)
        self.opened[relative][0].write(data)

    def close(self):
        for compressed, raw in self.opened.values():
            compressed.close(); raw.close()
        self.opened.clear()


def _preserve(review, destination, writer_limit):
    started, cpu = time.monotonic(), time.process_time()
    receipt = read_json(review)
    if receipt.get('schema') != REVIEW_SCHEMA:
        raise ValueError('AWARE preservation requires valid identification evidence')
    platforms = [p for p in ('illumina', 'ont', 'pacbio') if receipt.get('platforms', {}).get(p)]
    if not platforms:
        raise ValueError('No active input platform')
    for platform in set(('illumina', 'ont', 'pacbio')) - set(platforms):
        prefix = 'nanopore' if platform == 'ont' else platform
        for row in receipt.get('files', []):
            if Path(row['path']).name == prefix + '.kraken2.txt' and Path(row['path']).stat().st_size:
                raise ValueError('Inactive platform has Kraken assignments: ' + platform)
    from aware_bundle import LibraryWriter
    library_writer = LibraryWriter(destination, receipt.get('library_plan'))
    sample = receipt.get('isolate', '')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', sample):
        raise ValueError('Invalid isolate label')
    paths, before = {}, []
    required = {}
    if 'illumina' in platforms:
        required.update({k: v for k, v in FILES.items() if k.startswith('illumina') or k in ('r1', 'r2')})
    for platform in platforms:
        if platform != 'illumina':
            prefix = 'nanopore' if platform == 'ont' else 'pacbio'
            required.update({platform + '_report': prefix + '.report.txt', platform + '_calls': prefix + '.kraken2.txt', platform: platform + '_merged.fastq.gz'})
    for key, filename in required.items():
        rows = [r for r in receipt['files'] if Path(r['path']).name == filename]
        if len(rows) != 1:
            raise ValueError('Missing or duplicate review evidence: ' + filename)
        recorded = rows[0]; actual = identity(recorded['path'])
        if actual != recorded:
            raise ValueError('Reviewed input changed: ' + filename)
        paths[key] = Path(recorded['path']); before.append(actual)
    review_hash = digest(review)
    nodes, direct = taxonomy({p: paths[p + '_report'] for p in platforms})
    species_for, descendants = {}, {}
    for taxid in nodes:
        chain = lineage(nodes, taxid)
        species = next((x for x in chain if nodes[x]['rank'] == 'S'), None)
        if species:
            species_for[taxid] = species
            for ancestor in lineage(nodes, species):
                descendants.setdefault(ancestor, set()).add(species)
    bins = {}
    for taxid in sorted(set(species_for.values()), key=int):
        node = nodes[taxid]
        slug = re.sub(r'[^a-z0-9]+', '_', node['name'].lower()).strip('_')[:70] or 'taxon'
        bins[taxid] = dict(taxid=taxid, name=node['name'], rank='S',
            directory=sample + '_' + slug + '_taxid_' + taxid,
            category='provisional_species_assignment', units={p: 0 for p in platforms}, bases=0)
    for key in ('_ambiguous', '_unclassified'):
        bins[key] = dict(taxid=None, name=key[1:], directory=key, category=key[1:],
                         units={p: 0 for p in platforms}, bases=0)
    dbpath = destination / '.assignments.sqlite'
    db = sqlite3.connect(dbpath)
    db.execute('CREATE TABLE calls (platform TEXT, id TEXT, bucket TEXT, reason TEXT, candidates TEXT, unknown INTEGER, lengths TEXT, consumed INTEGER DEFAULT 0, PRIMARY KEY(platform,id))')
    observed = {p: Counter() for p in platforms}
    try:
        for platform in platforms:
            with paths[platform + '_calls'].open() as stream:
                for line in stream:
                    fields = line.rstrip('\r\n').split('\t')
                    if len(fields) != 5:
                        raise ValueError('Expected five-column Kraken assignments')
                    state, name, taxid, lengths, hits = fields
                    if state not in {'C', 'U'} or not taxid.isdigit() or (state == 'U') != (taxid == '0'):
                        raise ValueError('Invalid Kraken assignment status or TaxID')
                    expected = r'\d+\|\d+' if platform == 'illumina' else r'\d+'
                    if not re.fullmatch(expected, lengths) or any(int(n) < 1 for n in lengths.split('|')):
                        raise ValueError('Kraken read lengths do not match the platform')
                    if hits.split().count('|:|') != (1 if platform == 'illumina' else 0):
                        raise ValueError('Kraken paired-read hit separator does not match the platform')
                    bucket, reason, candidates, unknown = route(nodes, species_for, descendants, taxid, hits)
                    try:
                        db.execute('INSERT INTO calls VALUES (?,?,?,?,?,?,?,0)',
                            (platform, read_key(name, platform == 'illumina'), bucket, reason,
                             json.dumps(candidates), int(unknown), lengths))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError('Duplicate Kraken read identifiers; disambiguate the source reads') from exc
                    observed[platform][taxid] += 1
            if observed[platform] != direct[platform]:
                raise ValueError('Kraken assignment counts do not match the report: ' + platform)
        db.commit()
        writers = Writers(destination, writer_limit)
        reasons = Counter()
        try:
            with gzip.open(destination / 'assignments.tsv.gz', 'wt', newline='') as stream:
                audit = csv.writer(stream, delimiter='\t')
                audit.writerow(['platform', 'read_id', 'bucket', 'reason', 'candidate_species_taxids', 'unresolved_taxonomy'])
                for platform in platforms:
                    records = (itertools.zip_longest(fastq(paths['r1']), fastq(paths['r2']))
                               if platform == 'illumina' else ((r,) for r in fastq(paths[platform])))
                    for group in records:
                        if any(r is None for r in group):
                            raise ValueError('Illumina mate counts differ')
                        keys = [read_key(r[0], platform == 'illumina') for r in group]
                        if len(set(keys)) != 1:
                            raise ValueError('Illumina mate identifiers differ')
                        row = db.execute('SELECT bucket,reason,candidates,unknown,lengths,consumed FROM calls WHERE platform=? AND id=?', (platform, keys[0])).fetchone()
                        if row is None or row[-1]:
                            raise ValueError('Missing assignment or duplicate FASTQ identifier: ' + keys[0])
                        bucket, reason, candidates, unknown, lengths, _ = row
                        if list(map(int, lengths.split('|'))) != [r[2] for r in group]:
                            raise ValueError('FASTQ lengths differ from the Kraken assignment')
                        db.execute('UPDATE calls SET consumed=1 WHERE platform=? AND id=?', (platform, keys[0]))
                        info = bins[bucket]
                        library_writer.write(writers, info['directory'], platform, keys[0], group)
                        info['units'][platform] += 1; info['bases'] += sum(r[2] for r in group)
                        reasons[reason] += 1
                        audit.writerow([platform, keys[0], bucket, reason, ','.join(json.loads(candidates)), bool(unknown)])
            if db.execute('SELECT COUNT(*) FROM calls WHERE consumed=0').fetchone()[0]:
                raise ValueError('Kraken assignments contain reads absent from FASTQ inputs')
        finally:
            writers.close()
    finally:
        db.close()
    dbpath.unlink()
    for row in before:
        if identity(row['path'], full=False) != {k: v for k, v in row.items() if k != 'sha256'}:
            raise ValueError('Input changed during preservation')
    if digest(review) != review_hash:
        raise ValueError('Review receipt changed during preservation')
    for info in bins.values():
        info['linkage_status'] = 'not_evaluated' if sum(info['units'].get(p, 0) for p in ('ont', 'pacbio')) else 'no_assigned_long_reads'
        info['reconstruction_status'] = 'not_available'
        info['sample_manifest'] = library_writer.manifest(info['directory'], info['directory'])
    outputs = [{**identity(p), 'path': str(p.relative_to(destination))}
               for p in sorted(destination.rglob('*.gz'))]
    document = dict(schema=SCHEMA, version=VERSION, isolate=sample, created=datetime.now(timezone.utc).isoformat(),
        limitation=LIMITATION, review_sha256=review_hash, review_context=receipt.get('context'),
        source_sha256=digest(__file__), runtime={'python': runtime_platform.python_version(), 'system': runtime_platform.system()}, inputs=before, outputs=outputs,
        input_units={p: sum(c.values()) for p, c in observed.items()},
        output_units={p: sum(b['units'][p] for b in bins.values()) for p in observed},
        unit_definitions={p: 'read_pairs' if p == 'illumina' else 'reads' for p in platforms}, bins=list(bins.values()), reasons=dict(reasons),
        metrics={'wall_seconds': time.monotonic() - started, 'cpu_seconds': time.process_time() - cpu,
                 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if runtime_platform.system() == 'Darwin' else 1),
                 'scope': 'preservation process; peak RSS is the process lifetime high-water mark'})
    if document['input_units'] != document['output_units']:
        raise ValueError('Read conservation failed')
    (destination / 'manifest.json').write_text(json.dumps(document, indent=2) + '\n')
    rows = ''.join('<tr><td>' + html.escape(b['name']) + '</td><td>' + str(b['taxid'] or '—') +
                   '</td><td>' + str(b['units'].get('illumina', 0)) + '</td><td>' + str(sum(b['units'].get(p, 0) for p in ('ont', 'pacbio'))) + '</td></tr>'
                   for b in document['bins'])
    (destination / 'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1"><title>AWARE preservation</title>'
        '<style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:1rem;color:#18372b}th,td{padding:.6rem;border-bottom:1px solid #ccc;text-align:left}.scroll{overflow:auto}table{width:100%;border-collapse:collapse}</style>'
        '<h1>AWARE · Extended read preservation</h1><p>Experimental · ' + html.escape(sample) + '</p><p>' + LIMITATION + '</p>'
        '<p>Every input read or pair is retained in exactly one group. Shared candidates are recorded in '
        'assignments.tsv.gz; reads are not copied into multiple species groups. Empty groups have no FASTQ file.</p>'
        '<div class="scroll"><table><tr><th>Group</th><th>TaxID</th><th>Illumina pairs</th><th>Long reads</th></tr>' + rows +
        '</table></div><p><a href="manifest.json">Manifest and computational metrics</a> · <a href="assignments.tsv.gz">Assignment ledger</a></p></html>')
    return document


def preserve(review, output, writer_limit=32):
    """Publish an entire new output directory only after all integrity checks pass."""
    review, output = Path(review).absolute(), Path(output).absolute()
    if writer_limit < 1:
        raise ValueError('At least one writer is required')
    if output.exists() or output.is_symlink():
        raise ValueError('Output already exists; choose a new result directory')
    output.parent.mkdir(parents=True, exist_ok=True)
    lock = output.with_name(output.name + '.aware-plus.lock')
    with lock.open('x'):
        pass
    staging = None
    try:
        staging = Path(tempfile.mkdtemp(prefix='.' + output.name + '-', dir=output.parent))
        document = _preserve(review, staging, writer_limit)
        if output.exists() or output.is_symlink():
            raise ValueError('Output appeared during preservation; existing results are retained')
        staging.rename(output); staging = None
        return document
    finally:
        if staging is not None:
            shutil.rmtree(staging)
        lock.unlink()
