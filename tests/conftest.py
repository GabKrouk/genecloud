"""A tiny synthetic genome (60 genes) with the four source formats, so tests run offline in a second."""
import gzip
import random

import pytest

OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0006810
name: transport
namespace: biological_process
is_a: GO:0008150

[Term]
id: GO:0015706
name: nitrate transmembrane transport
namespace: biological_process
is_a: GO:0006810

[Term]
id: GO:0009408
name: response to heat
namespace: biological_process
is_a: GO:0008150
"""

NITRATE = "High affinity nitrate transporter involved in nitrate uptake from the soil."
HEAT = "Heat shock protein that acts as a molecular chaperone during heat stress."
FILLER = ["Pentatricopeptide repeat protein of the mitochondrion.", "Leucine-rich repeat receptor kinase.",
          "Zinc finger transcription factor.", "Cytochrome P450 monooxygenase.", "Glycosyl hydrolase.",
          "Ribosomal protein of the large subunit."]


def gid(i):
    return f"AT{1 + i % 5}G{10000 + i * 10:05d}"


@pytest.fixture(scope="session")
def sources(tmp_path_factory):
    d = tmp_path_factory.mktemp("src")
    rnd = random.Random(0)
    genes = [gid(i) for i in range(60)]
    text = {}
    for i, g in enumerate(genes):
        text[g] = NITRATE if i < 6 else HEAT if i < 10 else rnd.choice(FILLER)
    with gzip.open(d / "gene_info.gz", "wt") as fh:
        fh.write("#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\ttype_of_gene\t"
                 "Symbol_from_nomenclature_authority\tFull_name_from_nomenclature_authority\tNomenclature_status\t"
                 "Other_designations\tModification_date\tFeature_type\n")
        for i, g in enumerate(genes):
            sym = f"NRT{i}" if i < 6 else f"HSP{i}" if i < 10 else g
            fh.write(f"3702\t{i}\t{sym}\t{g}\t-\t-\t1\t-\t{text[g].split(' involved')[0]}\tprotein-coding\t-\t-\t-\t-\t20260101\t-\n")
    with gzip.open(d / "uniprot.tsv.gz", "wt") as fh:
        fh.write("Entry\tReviewed\tGene Names (ordered locus)\tGene Names (primary)\tGene Names (synonym)\tProtein names\t"
                 "Function [CC]\tKeywords\tProtein families\n")
        for i, g in enumerate(genes):
            kw = "Transport;Nitrate assimilation;Reference proteome" if i < 6 else "Stress response;Chaperone" if i < 10 else "Reference proteome"
            fh.write(f"P{i:05d}\treviewed\t{g.capitalize().replace('G', 'g')}\t\t\tProtein {i}\tFUNCTION: {text[g]} (PubMed:123) {{ECO:0000269}}.\t{kw}\t\n")
    with gzip.open(d / "annotations.gaf.gz", "wt") as fh:
        fh.write("!gaf-version: 2.2\n")
        for i, g in enumerate(genes):
            go = "GO:0015706" if i < 6 else "GO:0009408" if i < 10 else "GO:0008150"
            fh.write(f"TAIR\tlocus:{i}\tsym{i}\t\t{go}\tPMID:1\tIEA\t\tP\tname {i}\t{g}|sym{i}\tprotein\ttaxon:3702\t20260101\tTAIR\t\t\n")
    with gzip.open(d / "go-basic.obo.gz", "wt") as fh:
        fh.write(OBO)
    return d


@pytest.fixture(scope="session")
def annotation(sources, tmp_path_factory):
    from genecloud.annotation import build
    out = tmp_path_factory.mktemp("ann") / "toy.jsonl.gz"
    return build(sources, "arabidopsis", out)


@pytest.fixture(scope="session")
def engine(annotation):
    from genecloud import GeneCloud
    return GeneCloud(annotation, layers=("word", "phrase", "keyword", "go"), phrase_min_genes=3, cache=False)
