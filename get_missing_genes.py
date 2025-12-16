#!/usr/bin/env python3
import os
import sys
import pandas as pd

def die(msg: str, code: int = 1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def pick_first_glob(path_glob: str) -> str:
    import glob
    hits = glob.glob(path_glob)
    return hits[0] if hits else ""

def load_prokka_tsv(tsv_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    except Exception as e:
        die(f"Could not read TSV: {tsv_path} ({e})")

    required = {"locus_tag", "ftype", "length_bp", "gene", "product"}
    missing = required - set(df.columns)
    if missing:
        die(f"TSV missing expected columns {sorted(missing)}: {tsv_path}")

    # Keep CDS only
    df = df[df["ftype"].str.upper() == "CDS"].copy()

    # Normalize fields
    for c in ["locus_tag", "gene", "product"]:
        df[c] = df[c].astype(str).str.strip()

    # length_bp can be missing or non-numeric; keep as string but compute numeric where possible
    df["length_bp_num"] = pd.to_numeric(df["length_bp"], errors="coerce")

    return df

def is_hypothetical(prod: str) -> bool:
    p = (prod or "").lower()
    return "hypothetical protein" in p or p.strip() == ""

def main():
    ISO = sys.argv[1] if len(sys.argv) > 1 else "isolate1"

    # ---- Paths (adjust ONLY if your folders differ) ----
    hyb_tsv = pick_first_glob(f"{ISO}/annotation/*.tsv")
    spd_tsv = pick_first_glob(f"{ISO}/annotation_spades/*.tsv")

    if not hyb_tsv:
        die(f"Hybrid Prokka TSV not found at: {ISO}/annotation/*.tsv")
    if not spd_tsv:
        die(f"SPAdes Prokka TSV not found at: {ISO}/annotation_spades/*.tsv")

    out_dir = f"{ISO}/reports"
    os.makedirs(out_dir, exist_ok=True)

    out_csv = f"{out_dir}/missing_genes_recovery.csv"
    out_targets = f"{out_dir}/missing_targets.txt"

    print(f"--- S5: GENE RECOVERY AUDIT ({ISO}) ---")
    print(f"Hybrid TSV: {hyb_tsv}")
    print(f"SPAdes TSV: {spd_tsv}")

    hyb = load_prokka_tsv(hyb_tsv)
    spd = load_prokka_tsv(spd_tsv)

    print(f"Hybrid CDS rows:     {len(hyb)}")
    print(f"SPAdes CDS rows:     {len(spd)}")

    # Build SPAdes lookup sets
    spd_genes = set([g.lower() for g in spd["gene"].tolist() if g.strip()])
    spd_products = set([p.lower() for p in spd["product"].tolist() if p.strip()])

    missing_rows = []

    # Rule: for a civil engineer UX, we mimic the “what Prokka shows”
    # If gene symbol exists, compare by gene. Else compare by product (excluding hypothetical).
    for _, r in hyb.iterrows():
        gene = r["gene"].strip()
        prod = r["product"].strip()
        locus = r["locus_tag"].strip()
        length_bp = r["length_bp"].strip()

        missing = False
        missing_basis = ""

        if gene:
            if gene.lower() not in spd_genes:
                missing = True
                missing_basis = "missing_by_gene"
        else:
            if (not is_hypothetical(prod)) and (prod.lower() not in spd_products):
                missing = True
                missing_basis = "missing_by_product"

        if missing:
            missing_rows.append({
                "locus_tag_hybrid": locus,
                "gene_hybrid": gene if gene else "",
                "product_hybrid": prod if prod else "",
                "length_bp": length_bp,
                "basis": missing_basis
            })

    miss = pd.DataFrame(missing_rows)

    if miss.empty:
        print("No non-hypothetical CDS unique to HybAs found by Prokka report parity.")
        # Still write empty outputs for pipeline robustness
        miss.to_csv(out_csv, index=False)
        with open(out_targets, "w") as f:
            f.write("")
        print(f"Saved (empty): {out_csv}")
        print(f"Saved (empty): {out_targets}")
        return

    # Sort: longest first (usually more convincing in a figure/table)
    miss["length_bp_num"] = pd.to_numeric(miss["length_bp"], errors="coerce")
    miss = miss.sort_values(by=["length_bp_num"], ascending=False).drop(columns=["length_bp_num"])

    miss.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}  (rows={len(miss)})")

    # Targets for S6/S8:
    # Prefer locus_tag (most robust, always exists).
    # If locus_tag missing (rare), fall back to gene, then product.
    targets = []
    for _, r in miss.iterrows():
        lt = (r["locus_tag_hybrid"] or "").strip()
        gn = (r["gene_hybrid"] or "").strip()
        pr = (r["product_hybrid"] or "").strip()
        if lt:
            targets.append(lt)
        elif gn:
            targets.append(gn)
        elif pr:
            targets.append(pr)

    # unique, stable order
    seen = set()
    uniq = []
    for t in targets:
        if t not in seen:
            uniq.append(t)
            seen.add(t)

    with open(out_targets, "w") as f:
        for t in uniq:
            f.write(t + "\n")

    print(f"Saved: {out_targets}  (targets={len(uniq)})")
    print("Done.")

if __name__ == "__main__":
    main()
