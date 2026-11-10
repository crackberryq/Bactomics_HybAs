"""Presentation styles adapted from the BaSIS offline report."""
CSS = r"""
:root {
  color-scheme: light;
  --navy: #101f33;
  --navy-light: #1d334e;
  --ink: #1b2d40;
  --muted: #5e6c79;
  --paper: #f4f6f7;
  --white: #ffffff;
  --line: #dce3e8;
  --teal: #087c78;
  --teal-soft: #e8f5f2;
  --blue: #2858a2;
  --blue-soft: #eaf1fc;
  --amber: #89500b;
  --amber-soft: #fff3d8;
  --red: #a8303d;
  --red-soft: #fcecf0;
  --shadow: 0 3px 18px rgba(25, 47, 65, .04);
  --radius: 14px;
  font-synthesis: none;
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; scroll-padding-top: 24px; }
body { margin: 0; color: var(--ink); background: var(--paper); font: 15px/1.65 -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
a { color: var(--blue); text-underline-offset: 3px; }
a:hover { color: var(--teal); }
button, input, select { font: inherit; }
button, summary, select { -webkit-tap-highlight-color: transparent; }
button { cursor: pointer; }
:focus-visible { outline: 3px solid #d89925; outline-offset: 4px; border-radius: 4px; }
[hidden] { display: none !important; }
::selection { color: #102c36; background: #b7e5df; }
.skip-link { position: absolute; z-index: 100; top: -70px; left: 12px; color: #fff; background: var(--navy); padding: 10px 16px; border-radius: 6px; }
.skip-link:focus { top: 12px; }
.app { display: grid; grid-template-columns: 246px minmax(0, 1fr); min-height: 100vh; }
.sidebar { background: var(--navy); color: #d9e4ef; padding: 35px 23px 25px; height: 100vh; position: sticky; top: 0; display: flex; flex-direction: column; overflow-y: auto; }
.brand { display: flex; align-items: center; gap: 12px; color: #fff; text-decoration: none; margin: 0 0 35px; }
.brand:hover { color: #fff; }
.brand > svg { flex: 0 0 34px; }
.brand > span { display: block; font-size: 25px; font-weight: 800; letter-spacing: -.7px; line-height: 1.1; }
.brand > span > small { display: block; color: #acbece; font-size: 9px; font-weight: 650; letter-spacing: 1.1px; margin-top: 7px; }
.brand-mark { width: 40px; height: 40px; flex: 0 0 40px; display: grid; place-items: center; border: 1px solid #608b99; border-radius: 12px; color: #abe1d9; font-weight: 800; font-size: 22px; background: linear-gradient(140deg, #31596b, #173143); }
.brand-name { display: block; font-weight: 800; font-size: 25px; letter-spacing: -.7px; line-height: 1.15; }
.brand-subtitle { display: block; margin-top: 5px; color: #acbece; font-size: 10px; text-transform: uppercase; letter-spacing: 1.7px; font-weight: 650; }
.nav-label { margin: 0 11px 11px; color: #9aafc3; font-size: 10px; font-weight: 750; letter-spacing: 1.7px; text-transform: uppercase; }
.sidebar nav { display: flex; flex-direction: column; gap: 6px; }
.sidebar-caption { margin: 0 11px 12px; color: #9aafc3; font-size: 9px; font-weight: 750; letter-spacing: 1.25px; }
.nav-link { display: flex; align-items: center; gap: 12px; color: #bdcddd; text-decoration: none; font-size: 13px; padding: 12px; border: 1px solid transparent; border-radius: 9px; transition: background .15s, color .15s; }
.nav-link:hover, .nav-link.active, .nav-link[aria-current="location"] { color: #fff; background: #223b53; border-color: #34516a; }
.nav-link.active, .nav-link[aria-current="location"] { box-shadow: inset 3px 0 #69cdc1; }
.nav-link .nav-number { color: #93b6c8; font-size: 10px; font-family: ui-monospace, SFMono-Regular, Consolas, monospace; min-width: 19px; }
.nav-link > span { color: #93b6c8; font-size: 10px; font-family: ui-monospace, SFMono-Regular, Consolas, monospace; min-width: 19px; }
.sidebar-bottom { margin-top: auto; padding: 25px 10px 0; color: #d8e6ed; font-size: 11px; border-top: 1px solid #34495d; }
.sidebar-bottom small { display: block; color: #a8bbc9; font-size: 10px; line-height: 1.8; padding-top: 10px; }
.status-dot { display: inline-block; width: 6px; height: 6px; margin-right: 7px; border-radius: 50%; background: #71c7ad; }
.sidebar-foot { margin-top: auto; padding: 28px 10px 0; font-size: 11px; color: #a8bbc9; }
.sidebar-foot strong { display: block; color: #d8e6ed; font-size: 12px; font-weight: 650; margin-bottom: 6px; }
.sidebar-foot p { margin: 0 0 8px; }
.sidebar-rule { height: 1px; background: #34495d; margin-bottom: 19px; }
.main { min-width: 0; max-width: 1700px; padding: 0 42px 34px; width: 100%; }
.topbar { min-height: 81px; display: flex; justify-content: space-between; align-items: center; gap: 16px; color: var(--muted); font-size: 12px; }
.topbar .breadcrumb { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
.topbar strong { color: var(--ink); font-weight: 650; }
.topbar-actions { display: flex; align-items: center; justify-content: flex-end; gap: 9px; flex-wrap: wrap; }
.hero { display: flex; justify-content: space-between; align-items: center; gap: 24px; overflow: hidden; position: relative; isolation: isolate; color: #fff; padding: 40px 42px 35px; border-radius: 18px; background: radial-gradient(ellipse at 100% 0, #246276 0, transparent 48%), linear-gradient(115deg, #142d43, #173e52); }
.hero > div:first-child { min-width: 0; flex: 1 1 auto; }
.hero-art { flex: 0 0 250px; width: 250px; max-width: 32%; }
.hero-art svg { display: block; width: 100%; height: auto; }
.hero::after { content: ""; position: absolute; z-index: -1; height: 260px; width: 260px; right: -65px; top: -83px; border: 1px solid rgba(176, 231, 230, .17); border-radius: 50%; box-shadow: 0 0 0 32px rgba(176, 231, 230, .035), 0 0 0 65px rgba(176, 231, 230, .04), 0 0 0 99px rgba(176, 231, 230, .025); pointer-events: none; }
.eyebrow { margin: 0 0 10px; color: #39747b; font-size: 10px; font-weight: 750; letter-spacing: 2.1px; text-transform: uppercase; }
.hero .eyebrow { color: #a8dcd9; }
h1, h2, h3, h4 { color: inherit; font-weight: 730; line-height: 1.25; }
h1 { margin: 0 0 15px; font-size: clamp(28px, 2.8vw, 43px); letter-spacing: -1.25px; }
h2 { font-size: 24px; margin: 0; letter-spacing: -.6px; }
h3 { font-size: 16px; margin: 0 0 12px; letter-spacing: -.15px; }
h4 { font-size: 13px; margin: 20px 0 9px; }
p { margin: 0 0 13px; }
p:last-child { margin-bottom: 0; }
.hero-description, .hero > p:not(.eyebrow) { max-width: 830px; color: #dbe8ef; font-size: 14px; line-height: 1.8; }
.hero-meta { display: flex; flex-wrap: wrap; gap: 9px 24px; font-size: 12px; color: #c8dce5; margin-top: 24px; padding-top: 22px; border-top: 1px solid rgba(208, 231, 241, .2); }
.hero-meta strong { color: #fff; font-weight: 650; }
.hero .badge { vertical-align: middle; }
.section { margin-top: 36px; scroll-margin-top: 22px; }
.section-head { display: flex; align-items: flex-start; justify-content: flex-start; gap: 12px; margin-bottom: 18px; }
.section-head > div:first-child { min-width: 0; }
.section-head > .section-number { flex-shrink: 0; margin-top: 1px; margin-right: 0; }
.section-head > div { min-width: 0; }
.section-head .muted { font-size: 13px; margin-top: 7px; max-width: 940px; }
.section-title { display: flex; align-items: center; gap: 12px; }
.section-number { display: inline-flex; align-items: center; justify-content: center; min-width: 31px; height: 31px; border: 1px solid #bdccda; border-radius: 8px; color: #486982; font: 650 11px ui-monospace, SFMono-Regular, Consolas, monospace; margin-right: 8px; vertical-align: middle; }
.muted { color: var(--muted); }
.small { font-size: 12px; }
.mono, code { font-family: ui-monospace, SFMono-Regular, Consolas, "Liberation Mono", monospace; font-size: .86em; }
code { color: #36536b; background: #eef2f5; padding: 2px 5px; border-radius: 4px; overflow-wrap: anywhere; }
pre { overflow-x: auto; font-size: 11px; line-height: 1.65; background: #eef2f5; padding: 15px; border-radius: 8px; white-space: pre-wrap; overflow-wrap: anywhere; }
.stats { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 15px; margin: 20px 0; }
.stat { background: #fff; border: 1px solid var(--line); border-radius: var(--radius); padding: 23px 23px 20px; box-shadow: var(--shadow); min-width: 0; position: relative; overflow: hidden; }
.stat::before { content: ""; position: absolute; top: 0; left: 22px; width: 29px; height: 3px; background: var(--teal); border-radius: 0 0 4px 4px; }
.stat.warn::before { background: #b17613; }
.stat.info::before { background: var(--blue); }
.stat.neutral::before { background: #7c8c9b; }
.stat-label { color: #5a6c7b; font-size: 10px; line-height: 1.5; font-weight: 750; letter-spacing: 1.1px; text-transform: uppercase; }
.stat-value { display: block; margin: 10px 0 8px; color: var(--navy); font-size: clamp(27px, 2.7vw, 38px); line-height: 1.1; font-weight: 760; letter-spacing: -1.3px; font-variant-numeric: tabular-nums; overflow-wrap: anywhere; }
.stat-note { color: var(--muted); font-size: 11px; line-height: 1.6; }
.stat-value .unit { font-size: 15px; font-weight: 550; letter-spacing: 0; color: var(--muted); }
.grid-2 { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; align-items: start; }
.grid-3 { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 15px; }
.panel { min-width: 0; padding: 24px; background: #fff; border: 1px solid var(--line); border-radius: var(--radius); box-shadow: var(--shadow); margin-bottom: 18px; }
.panel-head { display: flex; align-items: flex-start; justify-content: space-between; gap: 12px; margin-bottom: 15px; }
.panel-head h3 { margin: 0; }
.panel .table-wrap { margin: 14px 0; border-width: 1px; border-radius: 9px; }
.panel > .table-wrap:only-child { margin: 0; }
.panel .table-wrap + p { margin-top: 13px; }
.table-wrap { overflow-x: auto; border: 1px solid var(--line); background: #fff; border-radius: var(--radius); }
table { width: 100%; border-collapse: collapse; font-size: 12px; text-align: left; }
caption { text-align: left; padding: 14px 20px; color: var(--muted); }
th { background: #f6f8fa; color: #5c6d7c; font-size: 10px; font-weight: 750; letter-spacing: .65px; text-transform: uppercase; padding: 13px 17px; border-bottom: 1px solid var(--line); vertical-align: bottom; white-space: nowrap; }
td { padding: 14px 17px; border-bottom: 1px solid #e9edf0; vertical-align: top; overflow-wrap: anywhere; }
tbody tr:last-child td { border-bottom: 0; }
tbody tr:hover { background: #f8fbfc; }
td .muted { display: block; margin-top: 3px; font-size: 11px; }
td code { font-size: 11px; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
.badge { display: inline-flex; align-items: center; width: fit-content; max-width: 100%; gap: 6px; padding: 4px 9px; color: #465a6c; background: #edf1f5; border: 1px solid #dce4ec; border-radius: 6px; font-size: 10px; font-weight: 720; line-height: 1.45; letter-spacing: .1px; overflow-wrap: anywhere; }
.badge.ok { color: #086854; background: #e8f5ef; border-color: #c0e4d4; }
.badge.warn { color: var(--amber); background: var(--amber-soft); border-color: #ecd79f; }
.badge.bad { color: #9a273b; background: var(--red-soft); border-color: #edc3ca; }
.badge.info { color: #275393; background: var(--blue-soft); border-color: #cbdaf3; }
.badge.neutral { color: #536372; background: #f0f3f5; border-color: #dbe2e8; }
.pill { display: inline-block; margin: 2px 4px 2px 0; padding: 3px 8px; background: #edf3f5; color: #4a6578; border: 1px solid #d9e3e8; border-radius: 5px; font-size: 10px; line-height: 1.5; }
.bar-track { height: 12px; border-radius: 6px; overflow: hidden; background: #e8edf1; display: flex; margin: 14px 0 12px; }
.bar-segment { min-width: 0; height: 100%; background: var(--teal); border-right: 2px solid #fff; }
.bar-segment:last-child { border-right: 0; }
.bar-segment.info { background: #376eb8; }
.bar-segment.warn { background: #c28b31; }
.bar-segment.neutral { background: #bac6d0; }
.bar-legend { display: flex; flex-wrap: wrap; gap: 8px 18px; color: var(--muted); font-size: 11px; }
.legend-dot { display: inline-block; width: 8px; height: 8px; border-radius: 2px; background: var(--teal); margin-right: 6px; }
.legend-dot.info { background: #376eb8; }
.legend-dot.warn { background: #c28b31; }
.legend-dot.neutral { background: #bac6d0; }
.score-display { display: flex; align-items: baseline; gap: 10px; flex-wrap: wrap; margin: 17px 0 9px; }
.score-value { color: var(--navy); font-size: 55px; line-height: 1.05; font-weight: 770; letter-spacing: -2.5px; font-variant-numeric: tabular-nums; }
.score-spread { font-size: 15px; color: var(--muted); }
.score-meter { height: 15px; background: linear-gradient(to right, #d5e9e5, #dce5ed); position: relative; border-radius: 5px; margin: 22px 0 7px; }
.score-marker { position: absolute; width: 4px; top: -5px; height: 25px; background: var(--navy); border: 1px solid white; border-radius: 3px; transform: translateX(-50%); }
.score-axis { display: flex; justify-content: space-between; color: var(--muted); font-size: 10px; font-variant-numeric: tabular-nums; }
.model-bars { display: grid; gap: 13px; margin: 16px 0; }
.model-row, .model-bar { display: grid; grid-template-columns: 65px minmax(0, 1fr) 60px; align-items: center; gap: 12px; color: var(--muted); font-size: 11px; }
.model-row .bar-track, .model-bar .bar-track { height: 8px; margin: 0; }
.model-row strong, .model-bar b { text-align: right; color: var(--ink); font-variant-numeric: tabular-nums; }
.model-result { display: flex; align-items: center; justify-content: space-between; gap: 25px; padding: 22px 24px; margin: 19px 0; background: linear-gradient(110deg, #eef7f5, #f4f7fa); border: 1px solid #d3e5e3; border-radius: 11px; }
.model-result > div:first-child { min-width: 0; }
.model-result p { font-size: 12px; }
.model-label { margin: 8px 0 10px; font-size: clamp(19px, 1.8vw, 25px); letter-spacing: -.5px; color: #245664; }
.model-score { flex: 0 0 auto; text-align: right; }
.model-score strong { display: block; color: #1a4356; font-size: clamp(35px, 3.3vw, 49px); font-weight: 760; line-height: 1.2; font-variant-numeric: tabular-nums; letter-spacing: -1.7px; }
.model-score span { display: block; font-size: 10px; color: #52717b; margin-top: 5px; }
.queue-summary { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 15px; padding: 8px 0 17px; margin-bottom: 13px; border-bottom: 1px solid #e1e8ed; }
.queue-summary > div { min-width: 0; }
.queue-summary span { display: block; font-size: 10px; color: var(--muted); }
.queue-summary strong { display: block; margin-top: 8px; color: #294f63; font-size: 28px; line-height: 1; letter-spacing: -.5px; }
.filters { display: flex; flex-wrap: wrap; gap: 10px; align-items: end; padding: 16px; margin-bottom: 16px; background: #eaf0f3; border: 1px solid #dce5eb; border-radius: 11px; }
.filters label { display: flex; flex-direction: column; gap: 5px; color: #4b6172; font-size: 10px; font-weight: 720; letter-spacing: .2px; }
.filters label:first-child { flex: 1 1 240px; }
.filters input, .filters select { min-height: 39px; min-width: 0; max-width: 100%; color: var(--ink); background: #fff; border: 1px solid #cbd7df; border-radius: 7px; padding: 8px 11px; font-size: 12px; font-weight: 450; letter-spacing: 0; }
.filters input { flex: 1 1 230px; }
.filters label input { flex: none; width: 100%; }
.filters select { cursor: pointer; }
.filters input::placeholder { color: #738594; opacity: 1; }
.filter-summary { display: flex; justify-content: space-between; align-items: center; gap: 14px; flex-wrap: wrap; margin: 13px 0; color: var(--muted); font-size: 12px; }
.filter-actions { display: flex; flex-wrap: wrap; gap: 8px; }
details.candidate { background: #fff; border: 1px solid var(--line); border-radius: 12px; margin-bottom: 10px; overflow: hidden; box-shadow: 0 2px 6px rgba(25, 47, 65, .02); }
details.candidate > summary { position: relative; display: flex; align-items: center; justify-content: space-between; gap: 18px; cursor: pointer; list-style: none; padding: 19px 48px 19px 22px; transition: background .15s; }
details.candidate > summary::-webkit-details-marker { display: none; }
details.candidate > summary::after { content: "+"; font-size: 22px; line-height: 1; font-weight: 450; position: absolute; top: 23px; right: 22px; color: #698096; }
details.candidate[open] > summary::after { content: "−"; }
details.candidate[open] > summary { background: #f4f8fa; border-bottom: 1px solid var(--line); }
details.candidate > summary:hover { background: #f4f8fa; }
.candidate-title { min-width: 0; font-size: 14px; font-weight: 700; line-height: 1.5; overflow-wrap: anywhere; }
.candidate-title .muted, .candidate-title small { display: block; font-size: 11px; font-weight: 450; margin-top: 3px; color: var(--muted); }
.candidate-title > span { display: block; font-size: 11px; font-weight: 450; margin-top: 3px; color: var(--muted); }
.candidate-meta { display: flex; flex-wrap: wrap; gap: 6px; justify-content: flex-end; align-items: center; font-size: 11px; min-width: 0; }
.candidate-body { padding: 24px; font-size: 12px; }
.candidate-body .grid-2 { gap: 23px; }
.candidate-body h3 { font-size: 14px; }
.candidate-body h4:first-child { margin-top: 0; }
.candidate-body .table-wrap { margin: 12px 0 20px; border-radius: 8px; }
.candidate-body table { font-size: 11px; }
.candidate-body th, .candidate-body td { padding: 10px 12px; }
.candidate-body .panel { padding: 18px; box-shadow: none; background: #fafcfd; }
.candidate-body ul { padding-left: 18px; }
.mini-stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(110px, 1fr)); gap: 10px; margin: 14px 0 20px; }
.mini-stats > div { background: #f5f8fa; border: 1px solid #e0e7ec; border-radius: 8px; padding: 11px 13px; color: var(--muted); font-size: 10px; }
.mini-stats strong { display: block; color: var(--ink); font-size: 18px; line-height: 1.4; font-weight: 700; font-variant-numeric: tabular-nums; }
.mini-stats .mini-label { display: block; color: var(--muted); font-size: 10px; }
.mini-stats > .stat { padding: 13px 15px; box-shadow: none; }
.mini-stats .stat::before { display: none; }
.mini-stats .stat-label { font-size: 9px; letter-spacing: .5px; }
.mini-stats .stat-value { font-size: 24px; margin: 6px 0 4px; }
.mini-stats .stat-note { font-size: 10px; }
.kv { display: block; font-size: 12px; margin: 0 0 16px; }
.kv > div { display: grid; grid-template-columns: minmax(100px, .85fr) minmax(0, 1.5fr); gap: 0 18px; border-bottom: 1px solid #e9edf1; }
.kv > div:last-child { border-bottom: 0; }
.kv dt, .kv dd { padding: 9px 0; margin: 0; border-bottom: 1px solid #e9edf1; overflow-wrap: anywhere; }
.kv > div > dt, .kv > div > dd { border-bottom: 0; }
.kv dt { color: var(--muted); }
.kv dd { font-weight: 550; }
.kv dt:last-of-type, .kv dd:last-of-type { border-bottom: 0; }
.evidence-list { list-style: none; margin: 0; padding: 0; }
.evidence-list li { position: relative; padding: 11px 0 11px 16px; border-bottom: 1px solid #e9edf1; overflow-wrap: anywhere; }
.evidence-list li::before { content: ""; position: absolute; left: 0; top: 20px; width: 5px; height: 5px; background: #87a6b2; border-radius: 50%; }
.evidence-list li:last-child { border-bottom: 0; }
details.evidence-list { padding: 0; margin: 15px 0; border: 1px solid #dce5eb; border-radius: 8px; background: #fafcfd; }
details.evidence-list > summary { cursor: pointer; padding: 12px 15px; color: #375c73; font-size: 11px; font-weight: 650; }
details.evidence-list[open] > summary { border-bottom: 1px solid #e1e8ee; margin-bottom: 12px; }
details.evidence-list > p, details.evidence-list > ul, details.evidence-list > dl, details.evidence-list > .table-wrap { margin-left: 15px; margin-right: 15px; margin-bottom: 15px; }
details.panel > summary { cursor: pointer; color: #36556b; font-size: 12px; }
details.panel[open] > summary { margin-bottom: 20px; }
.layer-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(185px, 1fr)); gap: 8px; margin: 12px 0 22px; }
.layer-grid > div { display: flex; align-items: center; justify-content: space-between; gap: 10px; min-width: 0; padding: 10px 12px; border: 1px solid #dce6eb; background: #f8fafb; border-radius: 7px; font-size: 10px; color: #4a6578; }
.layer-grid > div > span:first-child { flex: 1 1 auto; min-width: 0; }
.layer-grid .badge { font-size: 9px; flex-shrink: 0; }
.feature-map { display: block; width: 100%; max-width: 100%; height: auto; margin: 10px 0 14px; border: 1px solid #e1e8ee; border-radius: 8px; padding: 18px; background: #fafcfd; }
.feature-map text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; font-size: 11px; fill: #3e5c71; }
.feature-scroll { max-width: 100%; overflow-x: auto; margin: 10px 0 14px; border-radius: 8px; }
.feature-scroll .feature-map { min-width: 650px; margin: 0; }
.feature-actions { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; margin: 10px 0; }
.feature-actions p { margin: 0; color: #48677b; }
.feature-map [data-feature-detail] { cursor: pointer; }
.feature-map [data-feature-detail]:hover, .feature-map [data-feature-detail]:focus { stroke: #182b43; stroke-width: 2px; }
details.candidate:target { border-color: #228d83; box-shadow: 0 0 0 3px #cee8e2; }
.callout { padding: 16px 19px; margin: 16px 0; color: #39546c; background: #edf4f9; border: 1px solid #d4e2ed; border-left: 3px solid #6092b6; border-radius: 0 9px 9px 0; font-size: 12px; }
.callout strong { font-weight: 720; }
.callout.warn { color: #785119; background: #fff8e8; border-color: #e9d8b0; border-left-color: #c99738; }
.callout.ok { color: #276251; background: #edf7f2; border-color: #cfe4d9; border-left-color: #488f79; }
.callout.bad { color: #892e3d; background: #fff0f2; border-color: #efcbd3; border-left-color: #b9596a; }
.rec-card { margin-bottom: 13px; padding: 21px 23px; background: #fff; border: 1px solid var(--line); border-radius: 12px; box-shadow: var(--shadow); }
.rec-card > h3 { margin: 0 0 10px; }
.rec-card p { font-size: 12px; color: #4c6171; }
.rec-card .rec-meta { display: flex; align-items: center; flex-wrap: wrap; gap: 7px; margin-bottom: 13px; }
.rec-card .rec-trigger { font-size: 11px; background: #f3f6f8; border-radius: 7px; padding: 10px 12px; margin-top: 14px; }
.button { display: inline-flex; align-items: center; justify-content: center; gap: 7px; min-height: 35px; padding: 7px 12px; border: 1px solid #c8d5df; border-radius: 7px; background: #fff; color: #37546a; font-size: 11px; font-weight: 650; text-decoration: none; line-height: 1.5; transition: background .15s, border-color .15s; }
.button:hover { color: var(--navy); border-color: #8eabbc; background: #f1f6f8; }
.button-primary { color: #fff; background: var(--teal); border-color: var(--teal); }
.button-primary:hover { color: #fff; background: #086664; border-color: #086664; }
.empty { padding: 30px 24px; text-align: center; color: var(--muted); border: 1px dashed #c9d4dc; border-radius: 10px; background: #fff; font-size: 13px; }
.empty strong { color: var(--ink); display: block; margin-bottom: 7px; }
.audit-details { background: #f7f9fb; border: 1px solid #e0e6ed; border-radius: 8px; margin-top: 15px; }
.audit-details > summary { cursor: pointer; padding: 12px 15px; font-size: 12px; font-weight: 600; }
.audit-details > div { padding: 0 15px 15px; }
.download-list { display: flex; gap: 9px; flex-wrap: wrap; margin: 14px 0; }
.downloads { display: flex; gap: 9px; flex-wrap: wrap; margin: 14px 0; }
.text-link { font-size: 11px; font-weight: 650; white-space: nowrap; text-decoration: none; border-bottom: 1px solid #b1c6d9; padding-bottom: 2px; }
.path { overflow-wrap: anywhere; word-break: break-word; }
.status-line { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; margin: 11px 0; }
.section-footnote { color: var(--muted); font-size: 11px; margin: 13px 3px 0; }
.noscript { display: block; padding: 10px 14px; background: #fff4db; color: #7d561c; font-size: 12px; border: 1px solid #eed4a0; border-radius: 7px; margin-bottom: 14px; }
footer { border-top: 1px solid #d8e1e8; margin-top: 37px; padding: 23px 0 0; display: block; color: var(--muted); font-size: 11px; }
footer strong { color: #425b6e; font-weight: 650; }
@media (min-width: 1700px) { .main { padding-left: 56px; padding-right: 56px; } }
@media (max-width: 1250px) {
  .app { grid-template-columns: 214px minmax(0, 1fr); }
  .sidebar { padding-left: 17px; padding-right: 17px; }
  .main { padding-left: 26px; padding-right: 26px; }
  .hero { padding: 32px; }
  .hero-art { flex-basis: 210px; width: 210px; max-width: 29%; }
  .stats { gap: 11px; }
  .stat { padding: 21px 17px 18px; }
  .stat-value { font-size: 30px; }
  .grid-3 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
}
@media (max-width: 1020px) {
  .app { grid-template-columns: 188px minmax(0, 1fr); }
  .sidebar { padding-left: 12px; padding-right: 12px; }
  .brand { gap: 9px; margin-left: 4px; }
  .brand-name { font-size: 23px; }
  .brand-subtitle { font-size: 9px; letter-spacing: 1.1px; }
  .nav-link { font-size: 12px; padding: 11px 10px; gap: 8px; }
  .main { padding-left: 22px; padding-right: 22px; }
  .stats { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .grid-2 { grid-template-columns: minmax(0, 1fr); gap: 0; }
  .candidate-body .grid-2 { gap: 10px; }
  .hero { padding: 30px; }
  .hero-art { flex-basis: 150px; width: 150px; max-width: 24%; }
  .model-result { padding: 20px; gap: 18px; }
  .model-score strong { font-size: 37px; }
  .candidate-meta { max-width: 45%; }
}
@media (max-width: 740px) {
  html { scroll-padding-top: 14px; }
  .app { display: block; }
  .sidebar { height: auto; position: relative; padding: 19px 20px 12px; overflow: visible; }
  .brand { margin: 0 0 17px; }
  .brand-mark { width: 33px; height: 33px; flex-basis: 33px; font-size: 19px; border-radius: 9px; }
  .brand-name { font-size: 23px; }
  .brand-subtitle { font-size: 9px; }
  .sidebar nav { display: flex; flex-direction: row; overflow-x: auto; gap: 5px; padding-bottom: 5px; scrollbar-width: thin; scrollbar-color: #5b788e var(--navy); }
  .nav-link { flex: 0 0 auto; padding: 8px 10px; font-size: 11px; }
  .nav-number, .nav-label, .sidebar-foot, .sidebar-caption, .sidebar-bottom, .nav-link > span { display: none !important; }
  .main { padding: 0 17px 25px; }
  .topbar { min-height: 65px; font-size: 10px; align-items: center; }
  .topbar .breadcrumb { gap: 6px; }
  .topbar-actions { gap: 5px; }
  .hero { border-radius: 12px; padding: 26px 23px; }
  .hero-art, .desktop-only { display: none; }
  h1 { font-size: 29px; letter-spacing: -.7px; }
  h2 { font-size: 21px; }
  .hero-description, .hero > p:not(.eyebrow) { font-size: 12px; }
  .hero-meta { margin-top: 18px; padding-top: 17px; gap: 7px 15px; font-size: 10px; }
  .eyebrow { font-size: 9px; letter-spacing: 1.6px; }
  .section { margin-top: 29px; }
  .section-head { gap: 10px; flex-wrap: wrap; }
  .section-head .muted { font-size: 12px; }
  .section-number { min-width: 27px; height: 27px; font-size: 10px; }
  .stats { gap: 10px; margin: 15px 0; }
  .stat { padding: 19px 16px 17px; border-radius: 10px; }
  .stat-value { font-size: 29px; }
  .stat-label { font-size: 9px; letter-spacing: .65px; }
  .stat-note { font-size: 10px; }
  .panel { padding: 20px; border-radius: 11px; }
  .panel .table-wrap { margin: 13px 0; }
  .grid-3 { grid-template-columns: minmax(0, 1fr); }
  th, td { padding: 11px 13px; }
  .filters { padding: 12px; gap: 9px; }
  .filters label { flex: 1 1 130px; }
  .filters label:first-child { flex-basis: 100%; }
  .filters select { width: 100%; }
  details.candidate > summary { display: block; padding: 16px 39px 16px 16px; }
  details.candidate > summary::after { top: 17px; right: 16px; }
  .candidate-title { font-size: 13px; }
  .candidate-meta { max-width: none; justify-content: flex-start; margin-top: 10px; }
  .candidate-body { padding: 17px; }
  .mini-stats { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .mini-stats > div { padding: 9px 11px; }
  .kv { grid-template-columns: minmax(95px, .75fr) minmax(0, 1.4fr); gap: 0 11px; font-size: 11px; }
  .kv > div { grid-template-columns: minmax(95px, .75fr) minmax(0, 1.4fr); gap: 0 11px; }
  .model-result { display: block; padding: 18px; }
  .model-score { text-align: left; border-top: 1px solid #d3e3e5; padding-top: 13px; margin-top: 17px; }
  .model-score strong { font-size: 43px; }
  .model-label { font-size: 21px; }
  .layer-grid { grid-template-columns: minmax(0, 1fr); }
  .layer-grid > div { font-size: 11px; }
  .feature-map { padding: 8px; }
  .feature-map text { font-size: 15px; }
  .rec-card { padding: 19px; }
  .score-value { font-size: 46px; }
  footer { display: block; }
  footer > :last-child { margin-top: 8px; }
}
@media (prefers-reduced-motion: reduce) {
  html { scroll-behavior: auto; }
  *, *::before, *::after { transition: none !important; animation: none !important; }
}
.candidate-section-nav, .evidence-shortcuts { display:flex; flex-wrap:wrap; gap:8px 14px; }
.candidate-section-nav { margin-bottom:18px; }
.candidate-section { scroll-margin-top:24px; padding:12px 0; border-top:1px solid #dce6ec; }
.candidate-section:focus { outline:2px solid #147d7b; outline-offset:5px; }
.evidence-navigator td { vertical-align:top; }
.evidence-navigator .text-link { white-space:normal; }
.evidence-navigator td:nth-child(2) { min-width:170px; }
.evidence-navigator td:nth-child(4) { min-width:200px; }
@media print {
  @page { margin: 16mm 13mm; }
  html, body { background: #fff; color: #162b3d; font-size: 10px; }
  .app { display: block; }
  .sidebar, .topbar-actions, .filters, .filter-actions, .button, .skip-link, .noscript { display: none !important; }
  .main { padding: 0; max-width: none; }
  .topbar { min-height: 30px; }
  .hero { color: #162b3d; background: #f1f5f7 !important; border: 1px solid #c2d0db; padding: 22px; border-radius: 7px; break-inside: avoid; }
  .hero::after { display: none; }
  .hero-art { flex: 0 0 145px; width: 145px; max-width: 23%; }
  .hero p, .hero .hero-meta, .hero .hero-meta strong, .eyebrow { color: #345268 !important; }
  .hero-meta { border-color: #bfced8; }
  h1 { font-size: 27px; }
  h2 { font-size: 19px; }
  h3 { font-size: 13px; }
  h2, h3, h4, .section-head { break-after: avoid; }
  .section { margin-top: 22px; }
  .section-head { margin-bottom: 12px; }
  .stats { grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 9px; break-inside: avoid; }
  .stat { padding: 15px; box-shadow: none; }
  .stat-value { font-size: 25px; }
  .grid-2 { display: block; }
  .panel, .rec-card { box-shadow: none; padding: 17px; }
  .table-wrap { overflow: visible; }
  .feature-scroll { overflow: visible; }
  .feature-scroll .feature-map { min-width: 0; }
  table { font-size: 9px; }
  th { font-size: 8px; }
  th, td { padding: 8px 10px; }
  tr { break-inside: avoid; }
  thead { display: table-header-group; }
  .badge, .pill { font-size: 8px; }
  .callout, .rec-card, .stat { break-inside: avoid; }
  .bar-track, .bar-segment, .badge, .pill, .stat::before { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
  details.candidate { overflow: visible; box-shadow: none; margin-bottom: 12px; }
  details.candidate > summary { padding: 14px 17px; break-after: avoid; }
  details.candidate > summary::after { display: none; }
  .candidate-body { padding: 17px; }
  .candidate-title { font-size: 11px; }
  .candidate-meta { max-width: 50%; }
  a { color: inherit; text-decoration: underline; }
  pre { font-size: 8px; }
  footer { font-size: 9px; }
}
"""


CSS += r"""
.gx-search {padding:16px;background:var(--teal-soft);border-radius:10px;align-items:center}
#gx-overview {width:100%;cursor:crosshair;border:1px solid var(--line);border-radius:9px;background:#fafcfd}
.category-bar {display:grid;grid-template-columns:minmax(140px,1fr) 2fr 45px;gap:18px;align-items:center;width:100%;text-align:left;margin:8px 0;padding:12px 14px;text-transform:capitalize}
.category-track,.comparison-track {display:block;height:10px;background:#edf1f5;border-radius:5px;overflow:hidden;position:relative;min-width:90px}
.category-track>span {display:block;height:100%;background:var(--teal);border-radius:5px}
.comparison-track>span {position:absolute;height:100%;background:var(--teal);min-width:2px}
#gx-comparison td {vertical-align:middle}#gx-comparison small {display:block;font-size:10px}
#gx-comparison button {text-align:left}#gx-selection:empty {display:none}
@media(max-width:740px){.category-bar{grid-template-columns:minmax(90px,1fr) 1fr 30px;gap:8px}.gx-search{display:block}.gx-search input{margin:8px 0}}
"""
