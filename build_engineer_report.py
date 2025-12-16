#!/usr/bin/env python3
import os, sys, json, base64, re
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


def read_text(p: Path, default=""):
    try:
        return p.read_text(errors="replace")
    except Exception:
        return default


def b64img(p: Path):
    if not p or (not p.exists()):
        return None
    data = base64.b64encode(p.read_bytes()).decode("ascii")
    ext = p.suffix.lower().lstrip(".")
    mime = "image/png" if ext == "png" else "image/jpeg" if ext in ("jpg", "jpeg") else "application/octet-stream"
    return f"data:{mime};base64,{data}"


def load_config(cfg_path: Path):
    if cfg_path.exists():
        if cfg_path.suffix in (".yml", ".yaml"):
            if yaml is None:
                raise SystemExit("PyYAML not installed. Install with: pip install pyyaml")
            return yaml.safe_load(cfg_path.read_text())
        else:
            return json.loads(cfg_path.read_text())
    return {}


def html_escape(s: str):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;").replace("'", "&#039;"))


def find_first_existing(base: Path, candidates):
    for c in candidates:
        p = base / c if not isinstance(c, Path) else c
        if p.exists():
            return p
    return None


def find_figure(fig_dir: Path, stem: str):
    """Return the first matching figure with this stem (png/jpg/jpeg/svg)."""
    if not fig_dir.exists():
        return None
    exts = [".png", ".jpg", ".jpeg", ".svg"]
    for ext in exts:
        p = fig_dir / f"{stem}{ext}"
        if p.exists():
            return p
    # fallback: any file starting with stem.
    hits = list(fig_dir.glob(f"{stem}.*"))
    return hits[0] if hits else None


def parse_gff_attrs(attr: str):
    d = {}
    for part in attr.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def gff_lookup(gff_path: Path, locus_tag: str):
    """
    Find Name/gene/product for a locus_tag from a GFF.
    We match if:
      - 'locus_tag=LOCUS' in attributes OR
      - 'ID=LOCUS' in attributes
    """
    if not gff_path.exists():
        return ("", "", "")
    # scan linearly (GFF is not huge)
    with gff_path.open("r", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue
            attrs = parse_gff_attrs(parts[8])
            if attrs.get("locus_tag") == locus_tag or attrs.get("ID") == locus_tag:
                name = attrs.get("Name", "") or attrs.get("gene", "")
                gene = attrs.get("gene", "")
                product = attrs.get("product", "")
                return (name, gene, product)
    return ("", "", "")


def extract_verdict(proof_txt: Path):
    verdict = "INCONCLUSIVE"
    if not proof_txt.exists():
        return verdict
    for line in proof_txt.read_text(errors="replace").splitlines():
        if line.startswith("VERDICT:"):
            verdict = line.split("VERDICT:", 1)[1].strip()
            break
    return verdict


def main():
    if len(sys.argv) < 2:
        print("Usage: build_engineer_report.py <isolate> [final_root_or_audit_root] [config]")
        print("Tip: point it to FINAL_AUDIT_RESULTS so links (MultiQC/Evidence) work.")
        sys.exit(1)

    ISO = sys.argv[1]

    # If you pass FINAL_AUDIT_RESULTS, we build report inside FINAL_AUDIT_RESULTS/REPORT/
    # If you pass reports/full_audit, we build report inside reports/full_audit/
    ROOT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(ISO) / "reports" / "FINAL_AUDIT_RESULTS"
    CFG = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("engineer_report.yaml")

    cfg = load_config(CFG)
    title = cfg.get("title", "HybAs Engineer Audit Report")
    subtitle = cfg.get("subtitle", "Hybrid vs SPAdes: QC trajectory + annotation ghosting forensics")
    glossary = cfg.get("glossary", {})
    verdict_map = cfg.get("verdict_explanations", {})

    # Detect layout (FINAL bundle vs full_audit)
    # FINAL structure (recommended):
    #   FINAL_AUDIT_RESULTS/
    #     QC_ASSEMBLY/qc_steps
    #     FIGURES
    #     TABLES
    #     MULTIQC
    #     FORENSIC_GENES (optional flat copies)
    #     FORENSIC_EVIDENCE (full folders, recommended for clicking)
    #     GHOST_GENES (S5)
    #     REPORT (output here)
    is_final = (ROOT.name.upper() == "FINAL_AUDIT_RESULTS") or (ROOT / "REPORT").exists()

    if is_final:
        final_root = ROOT
        out_dir = final_root / "REPORT"
        out_dir.mkdir(parents=True, exist_ok=True)

        qc_dir = final_root / "QC_ASSEMBLY" / "qc_steps"
        tab_dir = final_root / "TABLES"
        fig_dir = final_root / "FIGURES"
        mq_dir = final_root / "MULTIQC"
        ghost_dir = final_root / "GHOST_GENES"
        evidence_dir = final_root / "FORENSIC_EVIDENCE"  # recommended to copy full folders here

        out_html = out_dir / "index.html"
        audit_label = str(final_root)
    else:
        # fallback: treat ROOT as audit root (reports/full_audit)
        audit_root = ROOT
        qc_dir = audit_root / "qc_steps"
        tab_dir = audit_root / "tables"
        fig_dir = audit_root / "figures"
        mq_dir  = audit_root / "multiqc"
        ghost_dir = Path(ISO) / "reports"
        evidence_dir = Path(ISO) / "forensic_evidence"

        out_html = audit_root / "index.html"
        audit_label = str(audit_root)

    # Tables
    benchmark_tsv = tab_dir / f"{ISO}_benchmark_table.tsv"
    reads_tsv     = tab_dir / f"{ISO}_reads_stats.tsv"

    # S5 outputs
    missing_csv = ghost_dir / "missing_genes_recovery.csv"
    targets_txt = ghost_dir / "missing_targets.txt"

    # Figures (auto-detect extension)
    uplot_fig = find_figure(fig_dir, f"{ISO}_benchmark_figure")
    adv_fig   = find_figure(fig_dir, f"{ISO}_advanced_metrics")
    gran_fig  = find_figure(fig_dir, f"{ISO}_granular_evolution")

    uplot_b64 = b64img(uplot_fig)
    adv_b64   = b64img(adv_fig)
    gran_b64  = b64img(gran_fig)

    # MultiQC html (pick any .html inside MULTIQC)
    mq_html = None
    if mq_dir.exists():
        # prefer isolate-named report if present
        preferred = list(mq_dir.glob(f"{ISO}*multiqc*.html"))
        if preferred:
            mq_html = preferred[0]
        else:
            htmls = list(mq_dir.glob("*.html"))
            mq_html = htmls[0] if htmls else None

    # Find GFF (to show gene name + product in ghost table)
    gff_path = None
    if not is_final:
        gffs = list((Path(ISO) / "annotation").glob("*.gff"))
        gff_path = gffs[0] if gffs else None
    else:
        # In FINAL folder, we *still* look back to isolate gff if available, but it's optional.
        gffs = list((Path(ISO) / "annotation").glob("*.gff"))
        gff_path = gffs[0] if gffs else None

    # Build ghost rows with gene name/product + verdict + evidence link
    ghost_rows = []
    if targets_txt.exists():
        targets = [x.strip() for x in targets_txt.read_text(errors="replace").splitlines() if x.strip()]
        for t in targets:
            # proof verdict: prefer FINAL evidence folder if exists, else original
            proof = (evidence_dir / t / "proof.txt") if evidence_dir else (Path(ISO) / "forensic_evidence" / t / "proof.txt")
            verdict = extract_verdict(proof)

            name, gene, product = ("", "", "")
            if gff_path:
                name, gene, product = gff_lookup(gff_path, t)

            explain = verdict_map.get(verdict, "See evidence bundle for details.")

            # Make evidence clickable (relative link)
            ev_folder = (evidence_dir / t) if evidence_dir else (Path(ISO) / "forensic_evidence" / t)
            ev_rel = os.path.relpath(ev_folder, out_html.parent) if ev_folder else ""

            ghost_rows.append({
                "locus": t,
                "name": name or gene or "",
                "product": product or "",
                "verdict": verdict,
                "explain": explain,
                "evidence_rel": ev_rel
            })

    # CSS / JS
    css = """
    body{font-family:Arial,Helvetica,sans-serif;margin:0;background:#0b0f14;color:#e6edf3}
    .wrap{max-width:1200px;margin:0 auto;padding:24px}
    .card{background:#111824;border:1px solid #243044;border-radius:14px;padding:16px;margin:14px 0}
    h1{margin:0 0 6px 0}
    .sub{opacity:.85;margin:0 0 10px 0}
    .tabs{display:flex;gap:8px;flex-wrap:wrap;margin:10px 0}
    .tabbtn{border:1px solid #2b3a52;background:#0f1622;color:#e6edf3;padding:10px 12px;border-radius:10px;cursor:pointer}
    .tabbtn.active{background:#1b2a40}
    .tab{display:none}
    .tab.active{display:block}
    table{width:100%;border-collapse:collapse}
    th,td{border-bottom:1px solid #243044;padding:8px;text-align:left;vertical-align:top}
    .small{font-size:13px;opacity:.88}
    img{max-width:100%;border-radius:12px;border:1px solid #243044}
    a{color:#7dd3fc}
    code{background:#0f1622;padding:2px 6px;border-radius:6px}
    .pill{display:inline-block;padding:2px 8px;border:1px solid #2b3a52;border-radius:999px;font-size:12px;opacity:.95}
    """

    js = """
    function showTab(id){
      document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
      document.querySelectorAll('.tabbtn').forEach(b=>b.classList.remove('active'));
      document.getElementById('tab_'+id).classList.add('active');
      document.getElementById('btn_'+id).classList.add('active');
    }
    """

    # Glossary
    gloss_html = ""
    if glossary:
        gloss_html += "<table><tr><th>Term</th><th>Meaning</th></tr>"
        for k, v in glossary.items():
            gloss_html += f"<tr><td><b>{html_escape(str(k))}</b></td><td>{html_escape(str(v))}</td></tr>"
        gloss_html += "</table>"
    else:
        gloss_html = "<p class='small'>No glossary provided in engineer_report.yaml</p>"

    # MultiQC block: iframe + fallback link
    mq_block = "<p class='small'>MultiQC file not found in MULTIQC/. (Install/run multiqc and re-run wrapper.)</p>"
    if mq_html and mq_html.exists():
        rel = os.path.relpath(mq_html, out_html.parent)
        mq_block = (
            f"<p class='small'>If the embedded view is blocked by your browser, open: "
            f"<a href='{html_escape(rel)}' target='_blank'>MultiQC report</a></p>"
            f"<iframe src='{html_escape(rel)}' style='width:100%;height:900px;border:1px solid #243044;border-radius:12px'></iframe>"
        )

    def img_block(label, b64, found_path: Path):
        if not b64:
            return f"<div class='card'><h3>{html_escape(label)}</h3><p class='small'>Missing figure in: <code>{html_escape(str(fig_dir))}</code></p></div>"
        note = html_escape(found_path.name) if found_path else ""
        return f"<div class='card'><h3>{html_escape(label)} <span class='pill'>{note}</span></h3><img src='{b64}'/></div>"

    plots_block = ""
    plots_block += img_block("S2: Benchmark trajectory (N50 + BUSCO)", uplot_b64, uplot_fig)
    plots_block += img_block("S3: Advanced dashboard (polish deltas + mapping + feature counts)", adv_b64, adv_fig)
    plots_block += img_block("S4: Granular evolution (BUSCO + N50 + deltas)", gran_b64, gran_fig)

    # Ghost table
    ghost_html = "<p class='small'>Targets produced by <code>S5</code> and validated by <code>S6→S8/S7</code>.</p>"
    if ghost_rows:
        ghost_html += (
            "<table><tr>"
            "<th>Target (locus_tag)</th>"
            "<th>Gene name</th>"
            "<th>Product</th>"
            "<th>Verdict</th>"
            "<th>Engineer explanation</th>"
            "<th>Evidence bundle</th>"
            "</tr>"
        )
        for r in ghost_rows:
            link = r["evidence_rel"]
            link_html = f"<a href='{html_escape(link)}'>open folder</a>" if link else "<span class='small'>not found</span>"
            ghost_html += (
                "<tr>"
                f"<td><b>{html_escape(r['locus'])}</b></td>"
                f"<td>{html_escape(r['name'])}</td>"
                f"<td class='small'>{html_escape(r['product'])}</td>"
                f"<td>{html_escape(r['verdict'])}</td>"
                f"<td>{html_escape(r['explain'])}</td>"
                f"<td><code>{html_escape(link)}</code><br/>{link_html}</td>"
                "</tr>"
            )
        ghost_html += "</table>"
    else:
        ghost_html += "<p>No targets found (missing_targets.txt missing/empty).</p>"

    # Tabs
    tab_def = [
        ("overview", "Overview", f"""
          <div class='card'>
            <h2>What this report answers</h2>
            <ul>
              <li>Did polishing improve structure (QUAST) and biology (BUSCO)?</li>
              <li>Which genes exist in HybAs but “ghost” in SPAdes annotation?</li>
              <li>Are ghost genes truly broken (frameshift) or just missed by annotation context?</li>
            </ul>
          </div>
          <div class='card'><h2>Glossary</h2>{gloss_html}</div>
        """),
        ("multiqc", "QC (MultiQC)", f"<div class='card'><h2>QC Summary</h2>{mq_block}</div>"),
        ("plots", "Plots", plots_block),
        ("ghosts", "Ghost genes table", f"<div class='card'><h2>Ghost Genes</h2>{ghost_html}</div>"),
    ]

    tab_buttons = ""
    tab_pages = ""
    for tid, name, body in tab_def:
        tab_buttons += f"<button class='tabbtn' id='btn_{tid}' onclick=\"showTab('{tid}')\">{html_escape(name)}</button>"
        tab_pages += f"<div class='tab' id='tab_{tid}'>{body}</div>"

    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{html_escape(title)} — {html_escape(ISO)}</title>
<style>{css}</style>
</head>
<body>
<div class="wrap">
  <div class="card">
    <h1>{html_escape(title)} — {html_escape(ISO)}</h1>
    <p class="sub">{html_escape(subtitle)}</p>
    <p class="small">Audit root: <code>{html_escape(audit_label)}</code></p>
  </div>

  <div class="tabs">{tab_buttons}</div>
  {tab_pages}

  <div class="card small">
    <b>Where outputs live (this bundle):</b><br/>
    <code>{html_escape(str(qc_dir))}</code> (QUAST/BUSCO per step)<br/>
    <code>{html_escape(str(tab_dir))}</code> (tables)<br/>
    <code>{html_escape(str(fig_dir))}</code> (figures)<br/>
    <code>{html_escape(str(mq_dir))}</code> (MultiQC)<br/>
    <code>{html_escape(str(evidence_dir))}</code> (per-gene evidence bundles)<br/>
  </div>
</div>

<script>{js}</script>
<script>showTab('overview');</script>
</body>
</html>
"""
    out_html.write_text(html)
    print(f"[+] Wrote engineer report: {out_html}")


if __name__ == "__main__":
    main()
