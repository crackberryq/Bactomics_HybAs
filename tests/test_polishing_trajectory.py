import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import hashlib

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("trajectory", ROOT/"workflow/scripts/polishing_trajectory.py")
trajectory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trajectory)


def test_retained_stage_report(tmp_path, monkeypatch):
    iso = tmp_path/"isolate"
    metadata = iso/"metadata"
    metadata.mkdir(parents=True)
    receipt = metadata/"polish_execution.json"
    stages = []
    for stage in ("draft","racon1","polypolish"):
        fasta = iso/(stage+".fa")
        fasta.write_text(">chr\nACGT\n")
        stages.append(dict(stage=stage,path=str(fasta),sha256=hashlib.sha256(fasta.read_bytes()).hexdigest()))
    receipt.write_text(json.dumps(dict(stages=stages, skipped_stages={"polypolish":"skipped"})))
    lineage = tmp_path/"lineage"
    lineage.mkdir()
    (lineage/"dataset.cfg").write_text("name=test")
    out = iso/"reports/polishing_trajectory"

    def run(cmd, check):
        if cmd[0] == "prodigal":
            faa = cmd[cmd.index("-a")+1]
            seq = "A" * (99 if faa.parent.name == "draft" else 100)
            faa.write_text(">a\n"+seq+"*\n")
            cmd[cmd.index("-o")+1].write_text("chr\tfixture\tCDS\t1\t4\t.\t+\t0\tID=a\n")
        elif cmd[0] == "busco":
            dest = cmd[cmd.index("--out_path")+1]
            table = dest/"busco/run_fixture/full_table.tsv"
            table.parent.mkdir(parents=True)
            table.write_text("a\tDuplicated\na\tDuplicated\nb\t"+("Fragmented" if dest.name=="draft" else "Complete")+"\n")
        elif cmd[0] == "dnadiff":
            prefix = cmd[2]
            prefix.with_suffix(".snps").write_text("2 C T 2 0 0 1 1 chr chr\n")
            prefix.with_suffix(".report").write_text("fixture")
        else:
            raise AssertionError(cmd)
    monkeypatch.setattr(trajectory.subprocess,"run",run)
    monkeypatch.setattr(trajectory.subprocess,"check_output",lambda *a,**kw:"fixture")
    monkeypatch.setattr(trajectory,"draw",lambda path,rows:path.write_text("<svg/>"))
    trajectory.measure(SimpleNamespace(receipt=receipt,out=out,lineage=lineage,threads=1))
    with (out/"transitions.tsv").open() as f:
        row, = csv.DictReader(f,delimiter="\t")
    assert row["response"]=="R" and float(row["delta_busco_pp"])==50
    assert float(row["relative_mean_change_pct"]) == 100/99
    assert row["aligned_substitutions"]=="1"
    assert "ID=a" in (out/"draft-to-racon1.edits.tsv").read_text()
    assert "b\tFragmented\tComplete" in (out/"marker_changes.tsv").read_text()
    assert not (out/"polypolish").exists()
    assert (out/"index.html").is_file()
