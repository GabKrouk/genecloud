"""Build a modern, per-gene annotation table from openly licensed sources.

Sources (see species.py): NCBI Gene (symbols, descriptions, nomenclature names), UniProtKB (protein names,
curated FUNCTION text, keywords, protein families), the Gene Ontology (GAF annotations + ontology) and, for
Arabidopsis, the latest TAIR public data release (Araport11 short descriptions, TAIR curator summaries,
computational descriptions, gene symbols and full names, mutant phenotypes, gene types).

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
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .ontology import Ontology
from .species import GO_OBO_URL, SPECIES, TAIR_API, TAIR_GENE_TYPE, TAIR_RELEASE_FILES, Species

csv.field_size_limit(sys.maxsize)

SOURCE_FILES = {
    "gaf": "annotations.gaf.gz",
    "obo": "go-basic.obo.gz",
    "gene_info": "gene_info.gz",
    "uniprot": "uniprot.tsv.gz",
}
TAIR_SOURCE_FILES = {
    "tair_descriptions": "tair_functional_descriptions.txt.gz",
    "tair_aliases": "tair_gene_aliases.txt.gz",
    "tair_phenotypes": "tair_germplasm_phenotypes.txt.gz",
    "tair_alleles": "tair_allele_phenotypes.txt.gz",
    "gene_type": "araport11_gene_type.txt.gz",
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
    tair_description: list = field(default_factory=list)    # Araport11 short descriptions (all gene models)
    tair_summary: str = ""                             # TAIR curator summary (literature-based, the richest text)
    tair_computational: list = field(default_factory=list)  # Araport11 computational descriptions
    phenotypes: list = field(default_factory=list)     # TAIR mutant phenotypes (germplasm / allele)
    biotype: str = ""                                  # Araport11 gene type (protein_coding, lncRNA, ...)


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
    if sp.tair:
        paths.update(download_tair(outdir, overwrite=overwrite))
    return paths


def _get(url: str, timeout: int = 600) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "genecloud/2"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def latest_tair_release() -> tuple[str, dict[str, str]]:
    """Find the most recent TAIR public data release. Returns (folder, {kind: file path})."""
    listing = json.loads(_get(f"{TAIR_API}/list?dir=Public_Data_Releases", 120))
    releases = sorted(k for k in listing if k.startswith("TAIR_Data_"))
    for rel in reversed(releases):
        files = {name: v["path"] for name, v in listing[rel].items() if isinstance(v, dict) and v.get("type") == "file"}
        found = {}
        for kind, pat in TAIR_RELEASE_FILES.items():
            hit = sorted(n for n in files if re.fullmatch(pat, n))
            if hit:
                found[kind] = files[hit[-1]]
        if "tair_descriptions" in found:
            return rel, found
    raise RuntimeError("no TAIR public data release with functional descriptions found")


def download_tair(outdir: str | Path, overwrite: bool = False) -> dict[str, Path]:
    """Download the latest TAIR public data release files (CC BY 4.0)."""
    outdir = Path(outdir)
    paths = {k: outdir / v for k, v in TAIR_SOURCE_FILES.items()}
    if all(p.exists() for p in paths.values()) and not overwrite:
        return paths
    rel, found = latest_tair_release()
    found["gene_type"] = TAIR_GENE_TYPE
    print(f"[genecloud] TAIR public data release: {rel}", file=sys.stderr)
    for kind, remote in found.items():
        p = paths[kind]
        if p.exists() and not overwrite:
            continue
        print(f"[genecloud] downloading {kind}: {remote}", file=sys.stderr)
        data = _get(f"{TAIR_API}/download?filePath={urllib.parse.quote(remote)}")
        if not data[:2] == b"\x1f\x8b":
            data = gzip.compress(data)
        p.write_bytes(data)
    (outdir / "tair_release.txt").write_text(rel + "\n")
    return {k: p for k, p in paths.items() if p.exists()}


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


_SOURCE_TAG = re.compile(r"\s*;?\s*\(source:[^)]*\)", re.I)
_NULL = {"", "NULL", "null", "-", "NA"}


def _tair_rows(path: Path):
    rd = csv.reader(_open_text(path), delimiter="\t")
    header = [h.strip().lower() for h in next(rd)]
    for row in rd:
        yield dict(zip(header, row))


def _locus(name: str, sp: Species, symbol_to_agi: dict[str, str]) -> str | None:
    m = sp.id_pattern.fullmatch(name.strip().split(".")[0]) if name else None
    if m:
        return sp.normalise(m.group(0))
    return symbol_to_agi.get(name.strip().upper())


def parse_tair_aliases(path: Path, sp: Species, genes: dict[str, Gene]) -> dict[str, str]:
    """Symbols and full names. Returns {SYMBOL: AGI} for symbols that point to a single locus."""
    owners: dict[str, set] = {}
    for r in _tair_rows(path):
        gid = _locus(r.get("locus_name", ""), sp, {})
        if not gid:
            continue
        g = genes.setdefault(gid, Gene(gid))
        sym, full = (r.get("symbol") or "").strip(), (r.get("full_name") or "").strip()
        if sym not in _NULL and not sp.id_pattern.fullmatch(sym):
            owners.setdefault(sym.upper(), set()).add(gid)
            if not g.symbol:
                g.symbol = sym
            elif sym != g.symbol and sym not in g.synonyms:
                g.synonyms.append(sym)
        if full not in _NULL and full not in g.names:
            g.names.append(full)
    return {s: next(iter(v)) for s, v in owners.items() if len(v) == 1}


def parse_tair_descriptions(path: Path, sp: Species, genes: dict[str, Gene]) -> None:
    """Araport11 functional descriptions, one row per gene model: merged per locus."""
    for r in _tair_rows(path):
        gid = _locus(r.get("name", ""), sp, {})
        if not gid:
            continue
        g = genes.setdefault(gid, Gene(gid))
        short = (r.get("short_description") or "").strip()
        if short not in _NULL and short not in g.tair_description:
            g.tair_description.append(short)
        summ = (r.get("curator_summary") or "").strip()
        if summ not in _NULL and summ not in g.tair_summary:
            g.tair_summary = (g.tair_summary + " " + summ).strip()
        comp = _SOURCE_TAG.sub("", (r.get("computational_description") or "")).strip(" ;")
        if comp not in _NULL and comp not in g.tair_computational:
            g.tair_computational.append(comp)
        gtype = (r.get("gene_model_type") or "").strip()
        if gtype not in _NULL and gtype != "unknown" and not g.biotype:
            g.biotype = gtype


def parse_tair_phenotypes(path: Path, sp: Species, genes: dict[str, Gene], symbol_to_agi: dict[str, str],
                          column: str = "phenotype", max_per_gene: int = 30) -> None:
    for r in _tair_rows(path):
        gid = _locus(r.get("locus_name", ""), sp, symbol_to_agi)
        if not gid or gid not in genes:
            continue
        ph = re.sub(r"\s+", " ", (r.get(column) or "")).strip()
        g = genes[gid]
        if ph not in _NULL and ph not in g.phenotypes and len(g.phenotypes) < max_per_gene:
            g.phenotypes.append(ph)


def parse_gene_type(path: Path, sp: Species, genes: dict[str, Gene]) -> None:
    for line in _open_text(path):
        if line.startswith("!") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        gid = _locus(f[0], sp, {})
        if gid and len(f) > 1 and gid in genes:
            genes[gid].biotype = f[1].strip()


def build(sources: dict[str, Path] | str | Path, species: str | Species = "arabidopsis",
          out: str | Path = "genecloud_annotation.jsonl.gz") -> Path:
    """Parse the downloaded sources and write the per-gene annotation table."""
    sp = SPECIES[species] if isinstance(species, str) else species
    if not isinstance(sources, dict):
        d = Path(sources)
        sources = {k: d / v for k, v in {**SOURCE_FILES, **TAIR_SOURCE_FILES}.items() if (d / v).exists()}
        rel = d / "tair_release.txt"
        tair_release = rel.read_text().strip() if rel.exists() else None
    else:
        tair_release = None
    onto = Ontology.from_obo(sources["obo"])
    genes: dict[str, Gene] = {}
    parse_gene_info(sources["gene_info"], sp, genes)
    parse_uniprot(sources["uniprot"], sp, genes)
    parse_gaf(sources["gaf"], sp, genes, onto)
    sym2agi: dict[str, str] = {}
    if "tair_aliases" in sources:
        sym2agi = parse_tair_aliases(sources["tair_aliases"], sp, genes)
    if "tair_descriptions" in sources:
        parse_tair_descriptions(sources["tair_descriptions"], sp, genes)
    if "tair_phenotypes" in sources:
        parse_tair_phenotypes(sources["tair_phenotypes"], sp, genes, sym2agi)
    if "tair_alleles" in sources:
        parse_tair_phenotypes(sources["tair_alleles"], sp, genes, sym2agi)
    if "gene_type" in sources:
        parse_gene_type(sources["gene_type"], sp, genes)
    meta = {"genecloud_annotation": 3, "species": sp.name, "taxon": sp.taxon, "built": _dt.date.today().isoformat(),
            "n_genes": len(genes), "sources": {k: str(Path(v).name) for k, v in sources.items()},
            "tair_release": tair_release}
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
            known = Gene.__dataclass_fields__
            genes[rec["id"]] = Gene(**{k: v for k, v in rec.items() if k in known})
    return genes, meta


def ontology_for(annotation_path: str | Path) -> Ontology | None:
    p = Path(annotation_path)
    obo = p.with_name(p.name.replace(".jsonl.gz", "") + ".go-basic.obo.gz")
    return Ontology.from_obo(obo) if obo.exists() else None
