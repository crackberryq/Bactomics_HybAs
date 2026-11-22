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


def test_density_preserves_order_and_distinguishes_unknown_rna():
    import base64
    import re
    from polishing_validation_figures import density_atlas
    windows=[{"contig":"c&d","start":"1","end":"10000","edit_bases":"1","rRNA_overlap_bp":""}]
    edits=[{"status":"mapped","event":"insertion","final_contig":"c&d","final_pos":"150","transition":"draft-to-final"}]
    fragment=density_atlas(windows,edits,[{"stage":"draft"},{"stage":"final"}])
    assert "c&amp;d" in fragment
    svg=base64.b64decode(re.search(r'base64,([^"]+)',fragment)[1]).decode()
    assert "Annotation coordinates unavailable" in svg
    assert "draft → final" in svg and "Insertion" in svg


def test_site_clusters_cover_every_locus_once():
    from polishing_validation import site_clusters
    assert list(site_clusters({("a",720),("a",199),("a",200),("b",201)}))==[
        ("a",[199,200]),("a",[720]),("b",[201])]


def test_render_only_rejects_stale_table_identity(tmp_path):
    import json
    import pytest
    from polishing_validation import render
    iso=tmp_path/"iso"
    traj=iso/"reports/polishing_trajectory"
    traj.mkdir(parents=True)
    tables=tmp_path/"tables"
    tables.mkdir()
    (traj/"stages.tsv").write_text("stage\tsha256\ndraft\t0000\n")
    (tables/"stage_to_final.tsv").write_text("stage\ndraft\n")
    (tables/"validation.json").write_text(json.dumps({"stages":{"draft":"1111"}}))
    with pytest.raises(ValueError,match="identities differ"):
        render(iso,tables,tmp_path/"edition")
    assert not (tmp_path/"edition/index.html").exists()


def test_complexity_dispersion_motifs_and_partial_windows():
    import pytest
    from polishing_validation_figures import spatial_context
    seq = {"r":"AAAAAACGTACGNN"}
    windows = [dict(contig="r",start=1,end=6,edit_bases=2,rRNA_overlap_bp="0"),
               dict(contig="r",start=7,end=12,edit_bases=0,rRNA_overlap_bp=""),
               dict(contig="r",start=13,end=14,edit_bases=1,rRNA_overlap_bp="1")]
    edits = [dict(transition="a-to-b",status="mapped",final_contig="r",final_pos=p) for p in [2,2,13]]
    support = [dict(contig="r",position=2,informative_bases=0),dict(contig="r",position=13,informative_bases=8),
               dict(contig="r",position=12,informative_bases=99)]
    complexity, dispersion, motifs, contexts = spatial_context(seq,windows,[windows[0],windows[2]],edits,support,["a-to-b"],6)
    assert [w["entropy_bits"] for w in complexity] == [0,2,""]
    assert dispersion[0]["variance_to_mean"] == 2
    assert dispersion[-1]["full_windows"] == 2
    assert motifs[0]["motif"] == "AAA" and motifs[0]["frequency_ratio"] > 1
    assert [c["annotation_context"] for c in contexts] == ["Outside annotated rRNA","Overlaps rRNA"]
    assert [c["mean_informative_depth"] for c in contexts] == [0,8]
    with pytest.raises(ValueError,match="Window boundaries"):
        spatial_context(seq,[dict(windows[0],edit_bases=9)],[],edits,[],["a-to-b"],6)


def test_refresh_replaces_inline_and_marked_legacy_sections_exactly_once():
    from polishing_validation import refresh_section
    old_inline = ('<h2>Polishing validation</h2><div>old inline values</div>'
                  '<p><a href="../polishing_validation/index.html">Download polishing validation tables</a></p>')
    old_backfill = ('<!-- HybAs final-coordinate validation --><section id="polishing-validation">'
                    '<h2>Polishing validation</h2>old backfill values'
                    '<a href="../polishing_validation/index.html">Validation evidence tables</a></section>')
    for inline, marked in [(old_inline,""),("",old_backfill),(old_inline,old_backfill)]:
        page = '<main><section id="stages">before'+inline+'</section><section>keep</section>'+marked+'</main>'
        result = refresh_section(page,"new values","validation/index.html")
        assert result.count("new values") == 1
        assert "old " not in result and "before" in result and "<section>keep</section>" in result
        updated = refresh_section(result,"newer values","validation2/index.html")
        assert updated.count("newer values") == 1 and "new values" not in updated
        assert refresh_section(updated,"newer values","validation2/index.html") == updated
