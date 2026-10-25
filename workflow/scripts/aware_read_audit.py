"""Disk-backed, read-only audit of AWARE read retention and Kraken assignments."""
import csv
import gzip
import os
import sqlite3
import tempfile
from pathlib import Path


def decision(path):
    if not path.is_file():
        return {}
    with path.open(encoding="utf8") as f:
        return {r["field"]: r["value"] for r in csv.DictReader(f, delimiter="\t") if r.get("field")}


def taxonomy(path):
    """Report indentation defines lineage. Clade counts are never added."""
    taxa, stack = {}, []
    if not path.is_file():
        return taxa
    with path.open(encoding="utf8", errors="replace") as f:
        for line in f:
            p = line.rstrip("\r\n").split("\t", 5)
            if len(p) != 6 or not p[4].strip().isdigit():
                continue
            depth = len(p[5]) - len(p[5].lstrip(" "))
            while stack and stack[-1][0] >= depth:
                stack.pop()
            lineage = [x[1] for x in stack] + [p[4].strip()]
            taxa[p[4].strip()] = {"name": p[5].strip(), "rank": p[3].strip(), "lineage": lineage}
            stack.append((depth, p[4].strip()))
    return taxa


def category(taxid, node, selected):
    if taxid == "0":
        return "Unclassified"
    if not node:
        return "Unresolved taxonomy"
    ids = set(node["lineage"])
    if "9606" in ids:
        return "Human-assigned"
    if "2" in ids:
        return "Selected bacterial lineage" if selected in ids else "Other bacteria"
    for tid, name in [("2157", "Archaea"), ("10239", "Viruses"), ("4751", "Fungi"),
                      ("2759", "Other eukaryotes")]:
        if tid in ids:
            return name
    return "Other classified taxa"


def canonical_id(header):
    word = header.strip().split()[0].removeprefix("@")
    return word[:-2] if word.endswith(("/1", "/2")) else word


def ids(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            if not h.startswith("@") or not plus.startswith("+") or not seq or not qual or len(seq.strip()) != len(qual.strip()):
                raise ValueError("Invalid FASTQ: " + str(path))
            yield canonical_id(h)


def audit(isolate, output_dir, selected_taxid):
    """Compare actual merged/cleaned read-pair IDs with saved Kraken assignments.

    This checks observed retention, not whether the classifier itself is correct.
    """
    before = isolate / "work/inputs_merged"
    after = isolate / "work/inputs_cleaned"
    kraken = isolate / "reports/kraken2"
    paths = {"b1": before / "illumina_R1.merged.fq.gz",
             "b2": before / "illumina_R2.merged.fq.gz",
             "a1": after / "illumina_R1.clean.fq.gz",
             "a2": after / "illumina_R2.clean.fq.gz"}
    assignments = kraken / "illumina.kraken2.txt"
    missing = [str(p) for p in [*paths.values(), assignments] if not p.is_file()]
    if missing:
        return {"status": "unavailable", "missing_files": missing}

    fd, dbfile = tempfile.mkstemp(prefix=".aware_audit_", suffix=".sqlite", dir=output_dir)
    os.close(fd)
    conn = sqlite3.connect(dbfile)
    try:
        conn.execute("CREATE TABLE reads(id TEXT PRIMARY KEY, b1 INT DEFAULT 0, b2 INT DEFAULT 0, "
                     "a1 INT DEFAULT 0, a2 INT DEFAULT 0, taxid TEXT, kc INT DEFAULT 0)")
        counts = {}
        for col, path in paths.items():
            sql = ("INSERT INTO reads(id," + col + ") VALUES (?,1) ON CONFLICT(id) DO UPDATE SET " +
                   col + "=" + col + "+1")
            batch, n = [], 0
            for rid in ids(path):
                batch.append((rid,))
                n += 1
                if len(batch) == 10000:
                    conn.executemany(sql, batch)
                    batch = []
            if batch:
                conn.executemany(sql, batch)
            conn.commit()
            counts[col] = n
        sql = ("INSERT INTO reads(id,taxid,kc) VALUES (?,?,1) ON CONFLICT(id) DO UPDATE SET "
               "kc=kc+1,taxid=CASE WHEN taxid IS NULL THEN excluded.taxid WHEN taxid=excluded.taxid THEN taxid ELSE '-1' END")
        batch, n = [], 0
        with assignments.open(encoding="utf8", errors="replace") as f:
            for line in f:
                p = line.rstrip("\r\n").split("\t")
                if len(p) < 3 or p[0] not in ("C", "U"):
                    continue
                taxid = p[2] if p[0] == "C" else "0"
                batch.append((canonical_id(p[1]), taxid if taxid.isdigit() else "-1"))
                n += 1
                if len(batch) == 10000:
                    conn.executemany(sql, batch)
                    batch = []
        if batch:
            conn.executemany(sql, batch)
        conn.commit()
        keys = ["input_pairs", "retained_pairs", "removed_pairs", "partial_pairs",
                "unpaired_input_ids", "unexpected_output_ids", "duplicate_fastq_ids",
                "unmatched_kraken_pairs", "duplicate_kraken_ids"]
        values = conn.execute("""
            SELECT
              SUM(b1=1 AND b2=1),
              SUM(b1=1 AND b2=1 AND a1=1 AND a2=1),
              SUM(b1=1 AND b2=1 AND a1=0 AND a2=0),
              SUM(b1=1 AND b2=1 AND a1!=a2),
              SUM(b1!=b2),
              SUM((a1>0 OR a2>0) AND (b1=0 OR b2=0)),
              SUM(b1>1 OR b2>1 OR a1>1 OR a2>1),
              SUM(b1=1 AND b2=1 AND kc=0),
              SUM(kc>1)
            FROM reads
        """).fetchone()
        totals = dict(zip(keys, [v or 0 for v in values]))
        taxa = taxonomy(kraken / "illumina.report.txt")
        by_taxon = []
        for taxid, b, r, removed, partial in conn.execute("""
            SELECT COALESCE(taxid,'-1'),COUNT(*),
              SUM(a1=1 AND a2=1),SUM(a1=0 AND a2=0),SUM(a1!=a2)
            FROM reads WHERE b1=1 AND b2=1 GROUP BY COALESCE(taxid,'-1')
        """):
            node = taxa.get(taxid)
            by_taxon.append({"taxid": taxid, "name": node["name"] if node else
                ("Unclassified" if taxid == "0" else "Unresolved assignment"),
                "rank": node["rank"] if node else "", "category": category(taxid, node, selected_taxid),
                "input_pairs": b, "retained_pairs": r or 0,
                "removed_pairs": removed or 0, "partial_pairs": partial or 0})
        by_taxon.sort(key=lambda row: (-row["removed_pairs"], row["taxid"]))
        categories = {}
        for row in by_taxon:
            bucket = categories.setdefault(row["category"], {"category": row["category"],
                "input_pairs": 0, "retained_pairs": 0, "removed_pairs": 0, "partial_pairs": 0})
            for k in ("input_pairs", "retained_pairs", "removed_pairs", "partial_pairs"):
                bucket[k] += row[k]
        valid = (totals["input_pairs"] == totals["retained_pairs"] + totals["removed_pairs"] and
                 all(totals[k] == 0 for k in keys[3:]))
        return {"status": "completed", "unit": "paired_read_ids", "input_file_records": counts,
                "kraken_assignment_rows": n, "totals": totals, "integrity_ok": valid,
                "by_category": sorted(categories.values(), key=lambda r: -r["removed_pairs"]),
                "by_taxon": by_taxon,
                "limitations": [
                    "Taxonomic matches are not independently verified biological origin.",
                    "Kraken paired-read assignments are reported per read-pair ID.",
                    "Missing report taxonomy is labelled unresolved; cumulative clade counts are not summed.",
                    "Observed retention is not proof that every whitelist/descendant choice was correct."
                ]}
    finally:
        conn.close()
        try:
            os.unlink(dbfile)
        except FileNotFoundError:
            pass


def audit_ont(isolate, output_dir, selected_taxid):
    """Single-read ONT counterpart; compare pre-AWARE merged and cleaned reads."""
    merged = isolate / "work/inputs_merged/ont_merged.fastq.gz"
    cleaned = isolate / "work/inputs_cleaned/ont_clean.fastq.gz"
    kr = isolate / "reports/kraken2/nanopore.kraken2.txt"
    report = isolate / "reports/kraken2/nanopore.report.txt"
    missing = [str(p) for p in (merged, cleaned, kr) if not p.is_file()]
    if missing:
        return {"status": "unavailable", "missing_files": missing}
    if merged.stat().st_size == 0:
        return {"status": "not_applicable", "reason": "No ONT input"}
    fd, dbfile = tempfile.mkstemp(prefix=".aware_ont_audit_", suffix=".sqlite", dir=output_dir)
    os.close(fd)
    conn = sqlite3.connect(dbfile)
    try:
        conn.execute("CREATE TABLE reads(id TEXT PRIMARY KEY, before_count INT DEFAULT 0, "
                     "after_count INT DEFAULT 0, taxid TEXT, kc INT DEFAULT 0)")
        for path, col in ((merged, "before_count"), (cleaned, "after_count")):
            sql = ("INSERT INTO reads(id," + col + ") VALUES (?,1) ON CONFLICT(id) DO UPDATE SET " +
                   col + "=" + col + "+1")
            batch = []
            for rid in ids(path):
                batch.append((rid,))
                if len(batch) >= 10000:
                    conn.executemany(sql, batch)
                    batch.clear()
            if batch:
                conn.executemany(sql, batch)
            conn.commit()
        batch = []
        sql = ("INSERT INTO reads(id,taxid,kc) VALUES (?,?,1) ON CONFLICT(id) DO UPDATE SET "
               "kc=kc+1,taxid=CASE WHEN taxid IS NULL THEN excluded.taxid WHEN taxid=excluded.taxid THEN taxid ELSE '-1' END")
        with kr.open(encoding="utf8", errors="replace") as f:
            for line in f:
                p = line.rstrip("\r\n").split("\t")
                if len(p) >= 3 and p[0] in ("C", "U"):
                    tid = p[2] if p[0] == "C" else "0"
                    batch.append((canonical_id(p[1]), tid if tid.isdigit() else "-1"))
                    if len(batch) >= 10000:
                        conn.executemany(sql, batch)
                        batch.clear()
        if batch:
            conn.executemany(sql, batch)
        conn.commit()
        names = ["input_reads", "retained_reads", "removed_reads", "unexpected_output_ids",
                 "duplicate_fastq_ids", "unmatched_kraken_reads", "duplicate_kraken_ids"]
        values = conn.execute("""
            SELECT SUM(before_count=1), SUM(before_count=1 AND after_count=1),
                   SUM(before_count=1 AND after_count=0),
                   SUM(after_count>0 AND before_count=0),
                   SUM(before_count>1 OR after_count>1),
                   SUM(before_count=1 AND kc=0), SUM(kc>1)
            FROM reads
        """).fetchone()
        totals = dict(zip(names, [v or 0 for v in values]))
        taxa = taxonomy(report)
        rows = []
        for tid, b, r, removed in conn.execute("""
            SELECT COALESCE(taxid,'-1'),COUNT(*),SUM(after_count=1),SUM(after_count=0)
            FROM reads WHERE before_count=1 GROUP BY COALESCE(taxid,'-1')
        """):
            node = taxa.get(tid)
            rows.append({"taxid": tid, "name": node["name"] if node else
                         ("Unclassified" if tid == "0" else "Unresolved assignment"),
                         "rank": node["rank"] if node else "",
                         "category": category(tid, node, selected_taxid),
                         "input_reads": b, "retained_reads": r or 0, "removed_reads": removed or 0})
        rows.sort(key=lambda x: (-x["removed_reads"], x["taxid"]))
        cats = {}
        for row in rows:
            c = cats.setdefault(row["category"], {"category": row["category"],
                               "input_reads": 0, "retained_reads": 0, "removed_reads": 0})
            for k in ("input_reads", "retained_reads", "removed_reads"):
                c[k] += row[k]
        valid = (totals["input_reads"] == totals["retained_reads"] + totals["removed_reads"]
                 and all(totals[k] == 0 for k in names[3:]))
        return {"status": "completed", "unit": "read_ids", "totals": totals,
                "integrity_ok": valid, "by_taxon": rows,
                "by_category": sorted(cats.values(), key=lambda x: -x["removed_reads"]),
                "limitations": ["Kraken assignments do not establish biological origin.",
                                "Taxonomic lineage is reconstructed from the saved report."]}
    finally:
        conn.close()
        try:
            os.unlink(dbfile)
        except FileNotFoundError:
            pass
