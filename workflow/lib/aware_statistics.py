"""Assignment-support uncertainty; never a calibrated probability of identity."""
from collections import Counter, defaultdict
import math
from pathlib import Path
import random
import tempfile

from aware_preservation import taxonomy, lineage, read_key
from aware_bundle import LibraryWriter


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return None
    if not 0 <= k <= n:
        raise ValueError('Invalid binomial counts')
    p = k / n
    denominator = 1 + z*z/n
    centre = (p + z*z/(2*n))/denominator
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/denominator
    return [max(0., centre-half), min(1., centre+half)]


def binomial(rng, n, p):
    """Exact beta-splitting sampler, with Bernoulli draws for small tails."""
    if p <= 0: return 0
    if p >= 1: return n
    offset = 0
    while n > 32:
        a = 1 + n//2
        b = n + 1 - a
        x = rng.betavariate(a, b)
        if x >= p:
            n, p = a-1, p/x
        else:
            offset += a
            n, p = b-1, (p-x)/(1-x)
    return offset + sum(rng.random() < p for _ in range(n))


def resample(rng, counts):
    remaining = sum(counts.values()); mass = remaining
    result = Counter()
    for key, count in sorted(counts.items()):
        value = binomial(rng, remaining, count/mass) if mass else 0
        result[key] = value
        remaining -= value; mass -= count
    return result


def write_report(path, nodes, counts):
    clades = Counter()
    for taxid, count in counts.items():
        if taxid != '0':
            for parent in lineage(nodes, taxid): clades[parent] += count
    n = sum(counts.values())
    rows = [f'{100*counts["0"]/n if n else 0}\t{counts["0"]}\t{counts["0"]}\tU\t0\tunclassified']
    for taxid, node in nodes.items():
        depth = len(lineage(nodes, taxid))-1
        rows.append(f'{100*clades[taxid]/n if n else 0:.10f}\t{clades[taxid]}\t{counts[taxid]}\t{node["rank"]}\t{taxid}\t'+ '  '*depth + node['name'])
    path.write_text('\n'.join(rows)+'\n' if n else '')


def analyse(receipt, selector, plan=None, replicates=200, seed=1729):
    paths = {Path(r['path']).name: Path(r['path']) for r in receipt['files']}
    platforms = [p for p in ('illumina','ont','pacbio') if receipt['platforms'].get(p)]
    prefixes = {p: 'nanopore' if p == 'ont' else p for p in platforms}
    nodes, totals = taxonomy({p: paths[prefixes[p]+'.report.txt'] for p in platforms})
    mapper = LibraryWriter('.', plan)
    strata = defaultdict(Counter)
    for p in platforms:
        seen = set()
        with paths[prefixes[p]+'.kraken2.txt'].open() as stream:
            for line in stream:
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) != 5:
                    raise ValueError('Expected five-column Kraken assignments for statistics')
                status, name, taxid, lengths, hits = fields
                name = read_key(name, p == 'illumina')
                if name in seen or taxid not in nodes and taxid != '0':
                    raise ValueError('Duplicate read or unresolved assignment taxonomy')
                seen.add(name)
                li, gi, lib, group = mapper.locate(p, name)
                # Independent resampling within each library/read group retains pairs.
                strata[(p, lib['library_id'], group['run_id'], group['read_group_id'])][taxid] += 1
        observed = Counter()
        for key, counts in strata.items():
            if key[0] == p: observed.update(counts)
        if observed != totals[p]:
            raise ValueError('Assignment counts disagree with reports; statistics withheld')
    target = receipt.get('recommendation', {}).get('taxid', '')
    def support(counts):
        n = sum(counts.values())
        k = sum(c for t,c in counts.items() if t != '0' and target in lineage(nodes,t))
        return dict(units=n, target_units=k, fraction=k/n if n else None, wilson_95=wilson(k,n))
    rows = [dict(platform=k[0], library=k[1], run=k[2], read_group=k[3], **support(v)) for k,v in sorted(strata.items())]
    by_platform = {p: support(totals[p]) for p in platforms}
    result = dict(target_taxid=target, target_scope='automatic AWARE recommendation before any manual override',
        unit='one Illumina read pair or one long read', libraries=rows, platforms=by_platform,
        interpretation='95% Wilson intervals describe observed assignment support under independent-unit sampling assumptions. Bootstrap selection frequency measures sensitivity to resampling these same assignments. Neither is a probability that the organism identification is correct. Database bias, duplicated molecules and other dependencies are not resolved; libraries and platforms are not pooled into a confidence score.')
    if selector and target:
        rng = random.Random(seed); votes = Counter()
        with tempfile.TemporaryDirectory(prefix='aware-statistics-') as directory:
            reports = {p: Path(directory)/(p+'.report.txt') for p in ('illumina','ont','pacbio')}
            for path in reports.values(): path.write_text('')
            for _ in range(replicates):
                aggregate = {p:Counter() for p in platforms}
                for key, counts in sorted(strata.items()): aggregate[key[0]].update(resample(rng,counts))
                for p in platforms: write_report(reports[p], nodes, aggregate[p])
                votes[selector(reports)['selected_taxid']] += 1
        result['stability'] = dict(method='stratified multinomial bootstrap of original Kraken assignment units',
            replicates=replicates, seed=seed, selections=dict(votes), same_target_fraction=votes[target]/replicates,
            policy='existing automatic target-election function, unchanged')
    return result
