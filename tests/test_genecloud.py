import re

from genecloud.text import Tokenizer, lemma
from tests.conftest import gid

NITRATE_SET = [gid(i) for i in range(6)]
BG = [gid(i) for i in range(60)]


def test_lemma():
    assert lemma("transporters") == "transporter"
    assert lemma("families") == "family"
    assert lemma("stress") == "stress"
    assert lemma("analysis") == "analysis"


def test_tokenizer_removes_noise():
    tok = Tokenizer(re.compile(r"\bAT[1-5]G\d{5}\b", re.I))
    words = tok.words("Nitrate transporters (PubMed:12345) {ECO:0000269} of AT1G01010 bind 3 ions; trans-Golgi.")
    assert "nitrate" in words and "transporter" in words
    assert not any(w.isdigit() or "pubmed" in w or "eco" in w or "at1g" in w for w in words)
    assert "trans‑golgi" in words                 # prefix compounds stay together


def test_build_and_load(annotation):
    from genecloud.annotation import load
    genes, meta = load(annotation)
    assert meta["n_genes"] == 60
    g = genes[gid(0)]
    assert "nitrate" in g.function.lower() and "PubMed" not in g.function
    assert "Transport" in g.keywords and g.go


def test_enrichment_finds_nitrate(engine):
    res = engine.run(NITRATE_SET, BG)
    top = res.significant
    assert len(top) >= 1
    sig = res.table[res.table.status == "significant"]
    assert any("nitrate" in c for c in sig.label)          # nitrate concepts may be merged under one representative
    row = res.table[res.table.label == "nitrate"].iloc[0]
    assert row.k == 6 and row.K == 6 and row.fdr < 0.01


def test_redundancy_merged(engine):
    res = engine.run(NITRATE_SET, BG)
    sig = res.table[res.table.status == "significant"]
    # all nitrate concepts are carried by the same 6 genes -> one representative
    nit = sig[sig.genes.str.count(",") == 5]
    assert nit.representative.sum() == 1


def test_min_genes_blocks_single_gene_words(engine):
    res = engine.run([gid(0), gid(20)], BG)
    assert (res.table[res.table.k < 2].fdr == 1).all()


def test_background_restricts_universe(engine):
    bg = [gid(i) for i in range(30)]
    res = engine.run(NITRATE_SET, background=bg)
    assert res.background_size == 30


def test_outputs(engine, tmp_path):
    res = engine.run(NITRATE_SET, BG)
    engine.plot(res, tmp_path / "c.png")
    engine.html(res, tmp_path / "c.html")
    res.to_tsv(tmp_path / "c.tsv")
    assert (tmp_path / "c.png").stat().st_size > 1000
    assert "<svg" in (tmp_path / "c.html").read_text()


def test_compare_and_legacy(engine, tmp_path):
    from genecloud.compare import draw_compare
    from genecloud.legacy import legacy_genecloud
    res = {"nitrate": engine.run(NITRATE_SET, BG), "heat": engine.run([gid(i) for i in range(6, 10)], BG)}
    d = draw_compare(res, tmp_path / "cmp.png")
    assert set(d.set) == {"nitrate", "heat"}
    t = legacy_genecloud(engine, NITRATE_SET, BG, sampling=20)
    assert t[t.word == "nitrate"].significant.iloc[0]


def test_cli(annotation, tmp_path):
    from genecloud.cli import main
    f = tmp_path / "genes.txt"
    f.write_text("# a comment line\n" + "\n".join(NITRATE_SET))
    b = tmp_path / "bg.txt"
    b.write_text("\n".join(BG))
    main(["run", str(f), "-a", str(annotation), "-b", str(b), "-o", str(tmp_path / "out"), "--formats", "png,html"])
    assert (tmp_path / "out.tsv").exists() and (tmp_path / "out.html").exists()


def test_background_is_mandatory(engine, annotation, tmp_path):
    import pytest
    from genecloud.cli import main
    with pytest.raises(TypeError):
        engine.run(NITRATE_SET)
    with pytest.raises(ValueError):
        engine.run(NITRATE_SET, None)
    f = tmp_path / "genes.txt"
    f.write_text("\n".join(NITRATE_SET))
    with pytest.raises(SystemExit):
        main(["run", str(f), "-a", str(annotation), "-o", str(tmp_path / "out")])
