"""Shared offline presentation for HybAs research evidence."""
from html import escape

STYLE = """
:root{--navy:#112d35;--ink:#1b2d40;--teal:#087c78;--paper:#f4f6f7;--line:#dce3e8;--muted:#5e6c79}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
.app{display:grid;grid-template-columns:230px minmax(0,1fr);min-height:100vh}
aside{background:var(--navy);color:#e1e9ec;padding:28px 18px;position:sticky;top:0;height:100vh}
.brand{font-size:27px;font-weight:750;color:white;margin-bottom:6px}.sub{font-size:10px;letter-spacing:1.5px;color:#9fc6c5;margin-bottom:32px}
aside a{display:block;color:#c7d8dc;text-decoration:none;padding:12px;border-radius:7px}aside a:hover{background:#25444d;color:white}
main{padding:30px 38px;max-width:1500px;width:100%}
.hero{background:linear-gradient(115deg,#142d43,#173e52);border-radius:17px;color:white;padding:35px;margin-bottom:24px}
.hero h1{font-size:34px;margin:0 0 8px}.hero p{margin:0;color:#c9e3e6}
.grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}
.stat,.panel{background:white;border:1px solid var(--line);border-radius:12px;padding:20px;margin-bottom:18px}
.stat small{color:var(--muted);display:block}.stat strong{font-size:25px;display:block;overflow-wrap:anywhere}
section{margin:34px 0}h2{font-size:23px}h3{font-size:16px}p.note{color:var(--muted)}
table{width:100%;border-collapse:collapse;font-size:12px}th{text-align:left;background:#edf3f2;color:#38565c}
td,th{border-bottom:1px solid var(--line);padding:10px;vertical-align:top;overflow-wrap:anywhere}
.table-wrap{overflow:auto;max-height:650px}tbody tr:hover{background:#f5f9f8}
code{overflow-wrap:anywhere}.badge{display:inline-block;background:#e6f2ef;color:#087c78;padding:4px 10px;border-radius:5px;font-size:12px}
.warn{background:#fff5df;color:#80540c;padding:12px;border-radius:8px}
input{width:100%;padding:10px;border:1px solid var(--line);border-radius:7px;margin-bottom:10px}
a{color:var(--teal)}footer{font-size:11px;color:var(--muted);padding:30px 0}
@media(max-width:850px){.app{display:block}aside{position:static;height:auto}main{padding:16px}.grid{grid-template-columns:repeat(2,1fr)}}

/* Bactologix-inspired refinement: clean hierarchy, responsive metrics and evidence graphs */
:root{--navy:#102c35;--teal:#136d68;--paper:#f5f7f6;--line:#dce4e3;--muted:#6d7e84}
body{font-size:14px;line-height:1.58}
aside{background:#112d35;padding:34px 19px}
aside a{font-size:13px;border-radius:9px;margin:4px 0}
main{padding:34px 44px;max-width:1520px}
.hero{background:radial-gradient(ellipse at 92% 0%,#276976 0%,transparent 52%),linear-gradient(115deg,#142d43,#173e52);padding:43px 40px;box-shadow:0 8px 28px #112d3510}
.hero h1{font-weight:720;letter-spacing:-1px}
.stat,.panel{box-shadow:0 3px 18px rgba(25,47,65,.045);border-radius:14px}
.stat{position:relative;overflow:hidden;padding:23px}
.stat:before{content:"";position:absolute;top:0;left:22px;width:30px;height:3px;background:#087c78}
.stat strong{font-size:clamp(22px,2.6vw,34px);font-weight:740;letter-spacing:-.8px}
section{margin-top:38px}section h2{font-weight:710;letter-spacing:-.5px}
.table-wrap{border:1px solid var(--line);border-radius:9px;margin-top:14px}
th{font-size:10px;text-transform:uppercase;letter-spacing:.7px;background:#f5f9f8}
td{font-variant-numeric:tabular-nums}
.stage-graph{padding:20px;background:#f8faf9;border:1px solid var(--line);border-radius:11px;margin:15px 0}
.stage-bar-row{display:grid;grid-template-columns:minmax(130px,1fr) minmax(100px,2fr) 85px;gap:13px;align-items:center;margin:11px 0;font-size:12px}
.stage-bar-row>span{overflow-wrap:anywhere}.stage-bar-row>strong{text-align:right;font-variant-numeric:tabular-nums}
.stage-track{height:15px;border-radius:5px;background:#e1e9e6;overflow:hidden}
.stage-fill{height:100%;border-radius:5px;background:linear-gradient(90deg,#137e79,#6bc1ad)}
@media(max-width:900px){main{padding:20px}.stage-bar-row{grid-template-columns:110px minmax(75px,1fr) 60px}}

.chart-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,410px),1fr));gap:12px}.chart-grid figure{margin:0}.chart-grid figcaption{padding-left:10px}.chart-grid svg{min-width:380px}
aside{overflow-y:auto}summary{cursor:pointer;font-weight:650}details{margin:12px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere}
.scroll{overflow:auto}main{min-width:0}img,svg{max-width:100%;height:auto}p,li{overflow-wrap:anywhere}
a:focus-visible,summary:focus-visible{outline:3px solid #14998e;outline-offset:3px}
.evidence-links{display:flex;flex-wrap:wrap;gap:9px 18px;padding:0;list-style:none}
.evidence-links a{display:inline-block;padding:5px 0}.stat small{font-size:12px}.panel h3:first-child{margin-top:0}
.evidence-chart{border:1px solid var(--line);border-radius:12px;overflow:auto;background:white;padding:16px;margin:14px 0}
.evidence-chart svg{display:block;width:100%;min-width:440px;height:auto}.panel>.table-wrap{margin-bottom:12px}
@media(max-width:850px){aside{height:auto;position:static;padding:16px 20px}.sub{margin-bottom:8px}aside nav{display:flex;flex-wrap:wrap;gap:2px}aside a{padding:5px 7px}.hero{padding:24px}.app{display:block}.grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media print{aside{display:none}.app{display:block}main{padding:0}.stat,.panel{box-shadow:none}.table-wrap,.scroll{overflow:visible;max-height:none}section{break-inside:auto}
.table-wrap table{font-size:9pt}.table-wrap th,.table-wrap td{padding:4pt;overflow-wrap:break-word}
.evidence-chart{overflow:visible;break-inside:avoid}.chart-grid{display:block}.chart-grid svg,.evidence-chart svg{min-width:0}
img,svg{break-inside:avoid}.panel img{max-height:210mm;object-fit:contain}}
"""

def evidence_page(title, context, body, nav=(), *, report_type="evidence"):
    marker = escape(report_type, quote=True)
    links = "".join(f'<a href="#{escape(anchor,quote=True)}">{escape(label)}</a>' for anchor,label in nav)
    return (f'<!doctype html><html lang="en" data-hybas-report="{marker}"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1"><title>HybAs · {escape(title)}</title>'
            f'<style>{STYLE}</style></head><body><div class="app">'
            '<aside><div class="brand">HybAs</div><div class="sub">BACTOMICS · RESEARCH EVIDENCE</div>'
            f'<nav aria-label="Report sections">{links}</nav></aside>'
            f'<main><header class="hero"><h1>{escape(title)}</h1><p>{escape(context)}</p></header>{body}'
            '<footer>Powered by GenoMix Human Link</footer></main></div></body></html>')
