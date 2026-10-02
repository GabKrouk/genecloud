"""Build a modern, per-gene annotation table from openly licensed sources.

Sources (see species.py): NCBI Gene (symbols, descriptions, nomenclature names), UniProtKB (protein names,
curated FUNCTION text, keywords, protein families) and the Gene Ontology (GAF annotations + ontology).

The result is one JSON record per gene, stored as gzipped JSON lines, with a small metadata header.
"""
from __future__ import annotations

import csv
import datetime as _dt
import gzip
import io
import json
import re
import sys
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .ontology import Ontology
from .species import GO_OBO_URL, SPECIES, Species

csv.field_size_limit(sys.maxsize)

SOURCE_FILES = {
    "gaf": "annotations.gaf.gz",
    "obo": "go-basic.obo.gz",
    "gene_info": "gene_info.gz",
    "uniprot": "uniprot.tsv.gz",
}


@dataclass
class Gene:
    id: str
    symbol: str = ""
    synonyms: list = field(default_factory=list)
    names: list = field(default_factory=list)          # full names (nomenclature, other designations)
    description: str = ""                              # NCBI Gene description
    protein_names: list = field(default_factory=list)  # UniProt recommended + alternative names
    function: str = ""                                 # UniProt FUNCTION comment (curated)
    keywords: list = field(default_factory=list)       # UniProt keywords
    families: list = field(default_factory=list)       # UniProt protein families
    go: list = field(default_factory=list)             # direct GO annotations (resolved ids)
    reviewed: bool = False                             # has a Swiss-Prot (reviewed) entry


# --------------------------------------------------------------------------------------------- download
def download(species: str | Species = "arabidopsis", outdir: str | Path = "genecloud_sources",
             overwrite: bool = False) -> dict[str, Path]:
    """Download the four source files for a species. Returns {kind: path}."""
    sp = SPECIES[species] if isinstance(species, str) else species
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    urls = {"gaf": sp.gaf_url, "obo": GO_OBO_URL, "gene_info": sp.gene_info_url, "uniprot": sp.uniprot_url}
    paths = {}
    for kind, url in urls.items():
        p = outdir / SOURCE_FILES[kind]
        paths[kind] = p
        if p.exists() and not overwrite:
            continue
        print(f"[genecloud] downloading {kind}: {url}", file=sys.stderr)
        req = urllib.request.Request(url, headers={"User-Agent": "genecloud/2"})
        with urllib.request.urlopen(req, timeout=600) as r:
            data = r.read()
        if not data[:2] == b"\x1f\x8b":                    # keep everything gzipped on disk
            data = gzip.compress(data)
        p.write_bytes(data)
    return paths


# ------------------------------------------------------------------------------------------------ parse
_PUBMED = re.compile(r"\s*\((?:PubMed|By similarity|Ref\.)[^)]*\)", re.I)
_ECO = re.compile(r"\s*\{ECO:[^}]*\}")
_PAREN_NAMES = re.compile(r"\(([^()]*(?:\([^()]*\))?[^()]*)\)")


def _open_text(path: Path):
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return io.StringIO(raw.decode("utf-8", errors="replace"))


def _split_protein_names(s: str) -> list[str]:
    """'Ammonium transporter 1 member 1 (AtAMT1;1) (Protein X)' -> ['Ammonium transporter 1 member 1', 'AtAMT1;1', 'Protein X']"""
    if not s:
        return []
    head = s.split(" (")[0].strip()
    rest = [m.strip() for m in _PAREN_NAMES.findall(s[len(head):])]
    rest = [r for r in rest if not r.startswith("EC ")]
    return [head] + rest


def _clean_function(s: str) -> str:
    s = s.replace("FUNCTION: ", " ")
    s = _ECO.sub("", s)
    s = _PUBMED.sub("", s)
    s = re.sub(r"\.(\s*\.)+", ".", s)
    return re.sub(r"\s+", " ", s).strip()


def parse_gene_info(path: Path, sp: Species, genes: dict[str, Gene]) -> None:
    rd = csv.reader(_open_text(path), delimiter="\t")
    header = next(rd)
    col = {c.lstrip("#"): i for i, c in enumerate(header)}
    for row in rd:
        m = sp.id_pattern.search(row[col["LocusTag"]]) or sp.id_pattern.search(row[col["Synonyms"]])
        if not m:
            continue
        gid = sp.normalise(m.group(0))
        g = genes.setdefault(gid, Gene(gid))
        sym = row[col["Symbol"]]
        if sym and sym != "-" and not sp.id_pattern.fullmatch(sym):
            g.symbol = g.symbol or sym
        syn = [x for x in row[col["Synonyms"]].split("|") if x not in ("-", "") and not sp.id_pattern.fullmatch(x)]
        g.synonyms += [x for x in syn if x not in g.synonyms]
        desc = row[col["description"]]
        if desc and desc != "-":
            g.description = desc
        for c in ("Full_name_from_nomenclature_authority", "Other_designations"):
            for x in row[col[c]].split("|"):
                if x and x != "-" and x not in g.names:
                    g.names.append(x)


def parse_uniprot(path: Path, sp: Species, genes: dict[str, Gene]) -> None:
    rd = csv.DictReader(_open_text(path), delimiter="\t")
    best: dict[str, list[dict]] = {}
    for row in rd:
        loci = {sp.normalise(m.group(0)) for m in sp.id_pattern.finditer(row.get("Gene Names (ordered locus)", "") or "")}
        for gid in loci:
            best.setdefault(gid, []).append(row)
    for gid, rows in best.items():
        rev = [r for r in rows if r.get("Reviewed") == "reviewed"]
        use = rev or rows
        g = genes.setdefault(gid, Gene(gid))
        g.reviewed = bool(rev)
        for r in use:
            for n in _split_protein_names(r.get("Protein names", "")):
                if n and n not in g.protein_names:
                    g.protein_names.append(n)
            f = _clean_function(r.get("Function [CC]", "") or "")
            if f and f not in g.function:
                g.function = (g.function + " " + f).strip()
            for k in (r.get("Keywords", "") or "").split(";"):
                k = k.strip()
                if k and k not in g.keywords:
                    g.keywords.append(k)
            fam = (r.get("Protein families", "") or "").strip()
            if fam and fam not in g.families:
                g.families.append(fam)
            prim = (r.get("Gene Names (primary)") or "").split()
            if prim and (not g.symbol or " " in g.symbol):      # NCBI symbols such as 'GSR 1' -> UniProt 'GLN1-1'
                g.symbol = prim[0]


def parse_gaf(path: Path, sp: Species, genes: dict[str, Gene], onto: Ontology) -> None:
    for line in _open_text(path):
        if not line or line[0] == "!":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or "NOT" in f[3]:
            continue
        m = sp.id_pattern.search(f[10]) or sp.id_pattern.search(f[2]) or sp.id_pattern.search(f[1])
        if not m:
            continue
        go = onto.resolve(f[4])
        if go is None:
            continue
        gid = sp.normalise(m.group(0))
        g = genes.setdefault(gid, Gene(gid))
        if go not in g.go:
            g.go.append(go)
        if f[9] and f[9] not in g.names and f[9] not in g.protein_names:
            g.names.append(f[9])


def build(sources: dict[str, Path] | str | Path, species: str | Species = "arabidopsis",
          out: str | Path = "genecloud_annotation.jsonl.gz") -> Path:
    """Parse the downloaded sources and write the per-gene annotation table."""
    sp = SPECIES[species] if isinstance(species, str) else species
    if not isinstance(sources, dict):
        d = Path(sources)
        sources = {k: d / v for k, v in SOURCE_FILES.items()}
    onto = Ontology.from_obo(sources["obo"])
    genes: dict[str, Gene] = {}
    parse_gene_info(sources["gene_info"], sp, genes)
    parse_uniprot(sources["uniprot"], sp, genes)
    parse_gaf(sources["gaf"], sp, genes, onto)
    meta = {"genecloud_annotation": 2, "species": sp.name, "taxon": sp.taxon, "built": _dt.date.today().isoformat(),
            "n_genes": len(genes), "sources": {k: str(Path(v).name) for k, v in sources.items()}}
    out = Path(out)
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        fh.write(json.dumps({"_meta": meta}) + "\n")
        for gid in sorted(genes):
            fh.write(json.dumps(asdict(genes[gid]), ensure_ascii=False) + "\n")
    # ship the ontology next to it (names are needed for the GO layer)
    obo_copy = out.with_name(out.name.replace(".jsonl.gz", "") + ".go-basic.obo.gz")
    obo_copy.write_bytes(Path(sources["obo"]).read_bytes() if str(sources["obo"]).endswith(".gz")
                         else gzip.compress(Path(sources["obo"]).read_bytes()))
    return out


def load(path: str | Path) -> tuple[dict[str, Gene], dict]:
    genes, meta = {}, {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            rec = json.loads(line)
            if "_meta" in rec:
                meta = rec["_meta"]
                continue
            genes[rec["id"]] = Gene(**rec)
    return genes, meta


def ontology_for(annotation_path: str | Path) -> Ontology | None:
    p = Path(annotation_path)
    obo = p.with_name(p.name.replace(".jsonl.gz", "") + ".go-basic.obo.gz")
    return Ontology.from_obo(obo) if obo.exists() else None
