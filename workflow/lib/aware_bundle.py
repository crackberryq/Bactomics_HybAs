"""Portable AWARE outputs, retained-library identity and explicit assembly handoff."""
from collections import defaultdict
from copy import deepcopy
import html
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import tempfile

from identification_review import digest, identity, read_json
from input_manifest import validate_manifest

SCHEMA = 'hybas-aware-bundle-v1'
CONTINUATION_SCHEMA = 'hybas-aware-continuation-v1'


class LibraryWriter:
    """Split namespaced merged inputs without losing library or chemistry metadata."""
    def __init__(self, root, plan=None):
        self.root, self.plan = Path(root), plan
        self.groups = defaultdict(dict)

    def locate(self, platform, name):
        if self.plan:
            libs = self.plan['libraries']
            if platform == 'illumina':
                libs = [x for x in libs if x['platform'] == platform]
                match = re.fullmatch(r'hybasL(\d+)G(\d+)R\d+', name)
            else:
                match = re.fullmatch(r'hybas_' + platform + r'_L(\d+)G(\d+)R\d+', name)
            if not match:
                raise ValueError('Cannot recover declared library for read: ' + name)
            li, gi = map(int, match.groups())
            try:
                lib = libs[li]; group = lib['read_groups'][gi]
            except IndexError as exc:
                raise ValueError('Read namespace is outside the declared library plan') from exc
            if lib['platform'] != platform:
                raise ValueError('Read namespace conflicts with its platform')
            return li, gi, lib, group
        if platform == 'pacbio':
            raise ValueError('PacBio preservation requires declared chemistry')
        lib = dict(library_id=platform, platform=platform)
        if platform == 'illumina':
            lib.update(library_type='paired_end', orientation='fr', preprocessing='fastp')
        return 0, 0, lib, dict(run_id='legacy', read_group_id='merged')

    def write(self, writers, bucket, platform, name, records):
        li, gi, lib, group = self.locate(platform, name)
        key = (platform, li, gi)
        entry = self.groups[bucket].setdefault(key, (deepcopy(lib), {k:v for k,v in group.items() if k not in {'r1','r2','reads'}}))
        keys = ('r1','r2') if platform == 'illumina' else ('reads',)
        for key, record in zip(keys, records):
            relative = f'{bucket}/{platform}_L{li}/G{gi}_{key}.fastq.gz'
            writers.write(relative, record[1])
            entry[1][key] = relative

    def manifest(self, bucket, sample_id):
        libraries = {}
        for (platform, li, gi), (lib, group) in self.groups[bucket].items():
            entry = libraries.setdefault((platform, li), {k:v for k,v in lib.items() if k != 'read_groups'})
            entry.setdefault('read_groups', []).append(group)
        if not libraries:
            return None
        if len(sample_id) > 100:
            sample_id = sample_id[:85] + '_' + hashlib.sha256(sample_id.encode()).hexdigest()[:12]
        sample = dict(sample_id=sample_id, libraries=list(libraries.values()), mode='auto')
        short = [l for l in sample['libraries'] if l['platform'] == 'illumina']
        if len(short) > 1 and len(short) == len(sample['libraries']):
            sample['mode'] = (self.plan or {}).get('selection', {}).get('selected_mode', (self.plan or {}).get('mode', 'dual_illumina_runs'))
            if sample['mode'] not in {'dual_illumina_runs', 'dual_illumina_complementary'}:
                sample['mode'] = 'dual_illumina_runs'
        return {'schema_version': 1, 'samples': [sample]}


def safe_file(root, relative):
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts or not p.parts:
        raise ValueError('AWARE bundle paths must be relative and stay within the bundle')
    root = Path(root).resolve()
    current = root
    for part in p.parts:
        current /= part
        if current.is_symlink():
            raise ValueError('AWARE bundle cannot contain symbolic links')
    if not current.is_file():
        raise ValueError('Missing AWARE bundle file: ' + relative)
    return current


def verify_bundle(path, full=True):
    path = Path(path)
    doc = read_json(path)
    if doc.get('schema') != SCHEMA:
        raise ValueError('Unsupported AWARE bundle')
    seen = set()
    for row in doc['files']:
        if row['path'] in seen:
            raise ValueError('Duplicate AWARE bundle file')
        seen.add(row['path'])
        p = safe_file(path.parent, row['path'])
        if p.stat().st_size != row['size'] or full and digest(p) != row['sha256']:
            raise ValueError('AWARE output changed: ' + row['path'])
    return doc


def verify_continuation(path, checksum, context, isolate):
    if not checksum or digest(path) != checksum:
        raise ValueError('Assembly continuation receipt changed')
    doc = read_json(path)
    if doc.get('schema') != CONTINUATION_SCHEMA:
        raise ValueError('Unsupported assembly continuation')
    bundle = Path(doc['bundle_path'])
    if digest(bundle) != doc['bundle_sha256']:
        raise ValueError('AWARE bundle changed after assembly was requested')
    data = verify_bundle(bundle)
    if data['context'] != context or data['isolate'] != isolate:
        raise ValueError('AWARE inputs, settings or source changed; create a new run')
    if not any(g['id'] == 'primary' and g.get('assembly_eligible') for g in data['groups']):
        raise ValueError('No supported primary read set is available for assembly')
    # Continuing uses the pipeline's existing cleaned inputs, so verify those too.
    for row in data['continuation_inputs']:
        if identity(row['path']) != row:
            raise ValueError('Saved AWARE inputs changed; cannot continue assembly')
    return data


def build(review, output, cleaned, decision, plan=None, preserve_groups=False, statistics=False, selector=None):
    from aware_preservation import Writers, fastq, read_key, preserve
    receipt = read_json(review)
    for expected in receipt['files']:
        if identity(expected['path']) != expected:
            raise ValueError('Identification evidence changed before AWARE completion')
    output = Path(output)
    if output.exists():
        raise ValueError('AWARE output already exists; retain it or create a new run')
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.aware-', dir=output.parent))
    try:
        mapper = LibraryWriter(staging, plan)
        writers = Writers(staging)
        primary_counts = {}
        before = [identity(p) for paths in cleaned.values() for p in paths] + [identity(decision)]
        try:
            for platform, paths in cleaned.items():
                count = 0
                for records in itertools.zip_longest(*(fastq(p) for p in paths)):
                    if any(r is None for r in records):
                        raise ValueError('Primary AWARE mate counts differ')
                    names = [read_key(r[0], platform == 'illumina') for r in records]
                    if len(set(names)) != 1:
                        raise ValueError('Primary AWARE mate identifiers differ')
                    mapper.write(writers, 'primary', platform, names[0], records)
                    count += 1
                primary_counts[platform] = count
        finally:
            writers.close()
        primary = mapper.manifest('primary', receipt['isolate'] + '_aware_primary')
        groups = [dict(id='primary', name='Primary AWARE output', category='primary', units=primary_counts, manifest=primary)]
        if preserve_groups:
            preservation = preserve(review, staging / 'groups')
            for group in preservation['bins']:
                manifest = group['sample_manifest']
                if manifest:
                    for lib in manifest['samples'][0]['libraries']:
                        for rg in lib['read_groups']:
                            for key in ('r1','r2','reads'):
                                if key in rg:
                                    rg[key] = 'groups/' + rg[key]
                groups.append(dict(id=group['directory'], name=group['name'], category=group['category'],
                                   taxid=group['taxid'], units=group['units'], manifest=manifest))
        for group in groups:
            group['assembly_eligible'] = False
            group['reason'] = 'No retained reads' if not group['manifest'] else 'Inspection/export only'
            if group['manifest'] and group['category'] in {'primary','provisional_species_assignment'}:
                try:
                    result = validate_manifest(group['manifest'], base_dir=staging)
                    group.update(assembly_eligible=True, mode=result['selections'][0]['selected_mode'], reason='Declared retained inputs are supported')
                except ValueError as exc:
                    group['reason'] = str(exc)
        stats = None
        if statistics:
            from aware_statistics import analyse
            stats = analyse(receipt, selector, plan)
            (staging/'statistics.json').write_text(json.dumps(stats, indent=2)+'\n')
        shutil.copyfile(decision, staging/'filtering.tsv')
        files = [dict(path=str(p.relative_to(staging)), size=p.stat().st_size, sha256=digest(p))
                 for p in sorted(staging.rglob('*')) if p.is_file()]
        doc = dict(schema=SCHEMA, isolate=receipt['isolate'], context=receipt['context'],
                   review_sha256=digest(review), parent_recommendation=receipt.get('recommendation'),
                   filtering_history='primary uses regular AWARE cleaning; species groups use provisional classification of original inputs',
                   continuation_inputs=before, files=files, groups=groups, statistics=stats,
                   limitation='Species groups are provisional and can contain multiple strains. Confidence intervals describe observed assignment support, not identification correctness.')
        (staging/'manifest.json').write_text(json.dumps(doc, indent=2)+'\n')
        esc=lambda x:html.escape(str(x))
        rows=''.join('<tr><td>'+esc(g['name'])+'</td><td>'+esc(g['category'])+'</td><td>'+esc(g['units'])+'</td><td>'+esc(g.get('mode', g['reason']))+'</td></tr>' for g in groups)
        (staging/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>AWARE outputs</title><style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:1rem}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #ccc}table{width:100%}</style><h1>HybAs AWARE · Outputs</h1><p>'+esc(doc['limitation'])+'</p><p>Select retained groups in Bactologix to register derived samples. Continue the parent run to assemble its primary output.</p><table><tr><th>Group</th><th>Category</th><th>Units by platform</th><th>Assembly route</th></tr>'+rows+'</table><p><a href="manifest.json">Portable bundle manifest</a> · <a href="filtering.tsv">Original filtering decisions</a></p>'+('<h2>Statistical evidence</h2><p>'+esc(stats['interpretation'])+'</p><pre>'+esc(json.dumps(stats,indent=2))+'</pre>' if stats else '')+'</html>')
        for row in before:
            if identity(row['path']) != row:
                raise ValueError('AWARE inputs changed during bundle creation')
        for expected in receipt['files']:
            observed = identity(expected['path'], full=False)
            if any(observed[k] != expected[k] for k in observed):
                raise ValueError('Identification evidence changed during AWARE completion')
        staging.rename(output)
        return doc
    finally:
        if staging.exists():
            shutil.rmtree(staging)
