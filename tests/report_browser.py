"""Desktop/mobile offline report navigation and export."""
import sys
import tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_report import fixture,write_report
from playwright.sync_api import sync_playwright
with tempfile.TemporaryDirectory() as tmp:
    root=Path(tmp);fixture(root);write_report(root,root/'reports/hybas')
    with sync_playwright() as pw:
        browser=pw.chromium.launch();page=browser.new_page(reduced_motion='reduce',viewport={'width':1440,'height':1000})
        errors=[];page.on('pageerror',lambda err:errors.append(str(err)))
        page.goto((root/'reports/hybas/report.html').as_uri())
        assert page.get_by_role('navigation').is_visible()
        page.get_by_role('link',name='AWARE decisions',exact=True).click()
        assert page.get_by_role('heading',name='Platform cleaning decisions').is_visible()
        assert page.get_by_role('cell',name='pass_through',exact=True).is_visible()
        assert page.get_by_role('cell',name='extract',exact=True).is_visible()
        assert page.locator('details[open]').count()==0
        page.get_by_role('link',name='Methods & provenance',exact=True).click()
        with page.expect_download() as d:page.get_by_role('button',name='Download report evidence').click()
        assert 'hybas-report-1' in Path(d.value.path()).read_text()
        with page.expect_download() as d:page.get_by_role('button',name='Download methods',exact=True).click()
        assert 'racon1, racon2, medaka' in Path(d.value.path()).read_text()
        page.get_by_role('link',name='Assembly overview',exact=True).click()
        Path('ui-previews').mkdir(exist_ok=True)
        page.screenshot(path='ui-previews/hybas-report-desktop.png',full_page=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
        page.screenshot(path='ui-previews/hybas-report-mobile.png',full_page=True)
        assert not errors,errors
        from test_aware_output_report import fixture as aware_fixture
        from aware_output_report import render
        aware=root/'aware.html';aware.write_text(render(aware_fixture()))
        page.set_viewport_size({'width':1440,'height':1000});page.goto(aware.as_uri())
        assert page.get_by_role('heading',name='Selected target and actual filtering').is_visible()
        assert page.get_by_role('cell',name='Manual choice · TaxID 22',exact=True).is_visible()
        page.get_by_role('link',name='Retained groups',exact=True).click()
        page.screenshot(path='ui-previews/aware-outputs-desktop.png',full_page=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
        page.screenshot(path='ui-previews/aware-outputs-mobile.png',full_page=True)
        assert not errors,errors
        browser.close()
