"""Durable identification evidence and explicit continuation; no policy changes."""
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path

SCHEMA = 'hybas-identification-review-v1'
DECISION_SCHEMA = 'hybas-identification-decision-v1'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def identity(path, full=True):
    p = Path(path)
    if p.is_symlink() or not p.is_file():
        raise ValueError('Identification evidence must be a regular file: ' + str(p))
    before = p.stat()
    result = dict(path=str(p.absolute()), size=before.st_size, mtime_ns=before.st_mtime_ns)
    if full:
        result['sha256'] = digest(p)
    after = p.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Identification evidence changed while reading: ' + str(p))
    return result


def context(config, sources, sample_sheet, inputs, databases):
    cfg = {k: v for k, v in config.items() if k not in {'kraken_review_decision', 'kraken_review_decision_sha256'}}
    doc = dict(config=cfg, sources={str(p): digest(p) for p in sources},
               sample_sheet=digest(sample_sheet) if sample_sheet else None,
               inputs=inputs, databases=databases)
    return hashlib.sha256(json.dumps(doc, sort_keys=True, default=str).encode()).hexdigest()


def read_json(path):
    p = Path(path)
    if p.is_symlink() or not p.is_file() or p.stat().st_size > 32 * 1024 * 1024:
        raise ValueError('Invalid or oversized identification review file: ' + str(p))
    return json.loads(p.read_text())


def load_decision(path, expected_hash):
    if not expected_hash or digest(path) != expected_hash:
        raise ValueError('Identification decision checksum differs from the recorded review')
    doc = read_json(path)
    if doc.get('schema') != DECISION_SCHEMA or doc.get('mode') not in {'aware', 'manual'}:
        raise ValueError('Unsupported identification decision')
    if not str(doc.get('taxid', '')).isdigit() or int(doc['taxid']) <= 1:
        raise ValueError('Select a valid identified taxon before continuing')
    return doc


def verify_decision(decision, current_context, isolate, full=True):
    path = decision['evidence_path']
    if digest(path) != decision['evidence_sha256']:
        raise ValueError('Identification report changed after review; review the evidence again')
    evidence = read_json(path)
    if evidence.get('schema') != SCHEMA or evidence.get('context') != current_context or evidence.get('isolate') != isolate:
        raise ValueError('Identification inputs, database, settings or source changed; create a new identification run')
    selected = next((x for x in evidence['candidates'] if x['taxid'] == decision['taxid']), None)
    if selected is None:
        raise ValueError('The selected TaxID is not in this identification report')
    if decision['mode'] == 'aware' and evidence['recommendation']['taxid'] != decision['taxid']:
        raise ValueError('The automatic decision does not match the recorded AWARE recommendation')
    for expected in evidence['files']:
        observed = identity(expected['path'], full=full)
        if any(observed[k] != expected[k] for k in observed):
            raise ValueError('Saved Kraken evidence changed: ' + expected['path'])
    return evidence


def write_report(path, html_path, isolate, context_id, rows, recommendation, files, include_children, platforms):
    doc = dict(schema=SCHEMA, isolate=isolate, created=datetime.now(timezone.utc).isoformat(),
               context=context_id, candidates=rows, recommendation=recommendation,
               include_children=include_children, platforms=platforms,
               files=[identity(p) for p in files])
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(doc, indent=2) + '\n')
    esc = lambda x: html.escape(str(x), quote=True)
    table = ''.join('<tr><td>' + esc(r['name']) + '</td><td>' + esc(r['taxid']) + '</td><td>' + esc(r['rank']) +
                    '</td><td>' + esc(r['illumina_pct']) + '</td><td>' + esc(r['ont_pct']) + '</td></tr>' for r in rows)
    Path(html_path).write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">'
        '<title>HybAs identification review</title><style>body{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:1rem;color:#18372b;background:#fff}'
        'table{border-collapse:collapse;width:100%}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #ccc}h1{color:#17664c}</style>'
        '<h1>HybAs AWARE · Identification review</h1><p>Sample: ' + esc(isolate) + '</p><h2>Awaiting researcher review</h2>'
        '<p>Filtering and assembly have not started. Return to Bactologix to choose a target and continue.</p>'
        '<p>AWARE recommendation: ' + esc(recommendation.get('name') or 'No supported candidate') + ' · TaxID ' + esc(recommendation.get('taxid') or '—') + '</p>'
        '<p>Percentages are Kraken clade assignments. Nested taxa overlap and must not be added. Illumina units are pairs; ONT units are reads. '
        'An absent platform is not evidence against a target. Classification is database-dependent.</p>'
        '<table><thead><tr><th>Taxon</th><th>TaxID</th><th>Rank</th><th>Illumina %</th><th>ONT %</th></tr></thead><tbody>' + table + '</tbody></table></html>')
    return doc
