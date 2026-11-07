"""Inspect the experimental inventory at desktop and phone widths."""
from pathlib import Path
import tempfile
from playwright.sync_api import sync_playwright
from test_aware_plus import evidence, preserve

with tempfile.TemporaryDirectory() as tmp:
    review, paths, output = evidence.__wrapped__(Path(tmp))
    preserve(review, output)
    with sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={'width': 1440, 'height': 1000})
        errors = []; page.on('pageerror', lambda err: errors.append(str(err)))
        page.goto((output/'report.html').as_uri())
        assert page.get_by_role('heading', name='AWARE+ · Read preservation').is_visible()
        assert page.get_by_role('cell', name='Example alpha', exact=True).is_visible()
        assert page.get_by_role('cell', name='ambiguous', exact=True).is_visible()
        assert 'not reconstructed isolates or guaranteed genomes' in page.locator('body').inner_text()
        assert page.get_by_role('link', name='Assignment ledger').get_attribute('href')=='assignments.tsv.gz'
        Path('ui-previews').mkdir(exist_ok=True)
        for width in (1440,390):
            page.set_viewport_size({'width':width,'height':900})
            assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
            page.screenshot(path=f'ui-previews/aware-plus-preservation-{width}.png',full_page=True)
        assert not errors,errors
        browser.close()
