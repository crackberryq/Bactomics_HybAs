"""Browser acceptance for the research report's saved polishing evidence."""
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from hybas_research_report import generate
from polishing_validation import render


def render_fixture(iso):
    trajectory = iso / 'reports/polishing_trajectory'
    validation = iso / 'reports/polishing_validation'
    trajectory.mkdir(parents=True)
    validation.mkdir(parents=True)
    final = iso / 'work/final/assembly.final.fasta'
    final.parent.mkdir(parents=True)
    final.write_text('>chr\n' + 'ACGT' * 5000 + '\n')
    hashes = {'draft': hashlib.sha256(b'fixture draft').hexdigest(),
              'polypolish': hashlib.sha256(final.read_bytes()).hexdigest()}
    (trajectory / 'stages.tsv').write_text(
        'stage\tsha256\tcomplete_pct\tshort_pct\tmean_aa\n' +
        f'draft\t{hashes["draft"]}\t96\t15\t250\npolypolish\t{hashes["polypolish"]}\t99\t10\t280\n')
    (validation / 'stage_to_final.tsv').write_text('stage\tsubstitutions\ndraft\t1\npolypolish\t0\n')
    (validation / 'validation.json').write_text(json.dumps(dict(stages=hashes, final_sha256=hashes['polypolish'],
        projected=1, unresolved=0, window_bp=10000, annotation='Fixture rRNA annotation', site_support='Fixture pileup')))
    header = 'contig\tstart\tend\tedit_bases\trRNA_overlap_bp\n'
    hotspot = 'chr\t1\t10000\t1\t10\n'
    (validation / 'windows.tsv').write_text(header + hotspot + 'chr\t10001\t20000\t0\t0\n')
    (validation / 'hotspots.tsv').write_text(header + hotspot)
    (validation / 'projected_edits.tsv').write_text(
        'status\tfinal_contig\tfinal_pos\tevent\ttransition\nmapped\tchr\t10\tsubstitution\tdraft-to-polypolish\n')
    (validation / 'site_support.tsv').write_text('contig\tposition\tinformative_bases\nchr\t10\t7\n')
    render(iso, validation, validation)
    research = iso / 'reports/hybas_research'
    generate(iso, research, False)
    return [validation / 'index.html', research / 'hybas_report.html']


def check_polishing(page, root):
    for report in render_fixture(root):
        for width, height in [(1440, 900), (390, 844)]:
            page.set_viewport_size({'width': width, 'height': height})
            page.goto(report.as_uri())
            atlas = page.locator('#edit-density-atlas img')
            assert atlas.count() == 1 and atlas.is_visible()
            assert atlas.evaluate('el => el.complete && el.naturalWidth > 0')
            page.get_by_text('Polishing trajectory projections', exact=True).click()
            projection = page.get_by_alt_text('Top, front and side projections', exact=False)
            assert projection.is_visible() and projection.evaluate('el => el.naturalWidth > 0')
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
            assert page.locator('nav a[href^="#"], aside a[href^="#"]').evaluate_all(
                'links => links.length > 0 && links.every(a => document.getElementById(a.hash.slice(1)))')
            atlas.scroll_into_view_if_needed()
            if report.name == 'index.html':
                page.screenshot(path=f'ui-previews/polishing-density-{width}.png')
