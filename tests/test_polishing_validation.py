import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"workflow/scripts"))
from polishing_validation import annotation, overlap_bp, projection


def test_reverse_gap_and_ambiguous_projection(tmp_path):
    path=tmp_path/"test.delta"
    path.write_text("ref query\nNUCMER\n>r q 8 8\n1 4 5 1 0 0 0\n-3\n0\n>r2 q 8 8\n1 1 1 1 0 0 0\n0\n")
    hits=projection(path,{"q":[1,2,3,4,5]})
    assert hits=={("q",1):None,("q",2):("r",3),("q",4):("r",2),("q",5):("r",1)}


def test_reference_provenance_and_overlap(tmp_path):
    gff=tmp_path/"a.gff"
    gff.write_text("other\tt\trRNA\t2\t4\t.\t+\t.\tID=r\n##FASTA\n>other\nACGT\n")
    counts,features,status=annotation(gff,{"ctg":"ACGT"})
    assert counts["rRNA"]==1 and features["ctg"]==[(2,4)]
    assert "verified" in status
    assert "unresolved" in annotation(gff,{"ctg":"ACTT"})[2]
    assert overlap_bp([(2,4),(3,6)],3,5)==3
