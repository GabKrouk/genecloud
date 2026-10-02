"""Export a GeneCloud annotation as compact static files for the web application (``web/``).

The web app runs the whole analysis in the browser (no server): it downloads these files once, then the
hypergeometric test, Benjamini-Hochberg correction, redundancy merging and the cloud layout are computed in
JavaScript (``web/engine.js``), exactly as in :meth:`genecloud.core.GeneCloud.run`.

Files written to ``outdir``:
  meta.json          annotation metadata, layers, text presets, background presets, file list
  genes.json.gz      gene ids, symbols, short descriptions, flags (ATH1 array, protein coding), per-layer
                     "annotated" bits
  concepts_<x>.json.gz   labels of the concepts of one block
  index_<x>.bin.gz   members of each concept of one block: varint(count), then delta-encoded gene indices

Blocks: ``wp_all``, ``wp_tair``, ``wp_open`` (words + phrases for each text preset), ``kw`` (UniProt keywords)
and ``go`` (GO terms with their ancestors). Concepts carried by fewer than two genes of the whole annotation
cannot be tested (min_genes >= 2) and are not exported; the per-gene layer bits keep the background exact.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

from .core import TEXT_PRESETS, GeneCloud

LAYER_CODE = {"word": "w", "phrase": "p", "keyword": "k", "GO BP": "b", "GO MF": "f", "GO CC": "c", "GO ": "g"}
LAYER_BIT = {"w": 1, "p": 2, "k": 4, "b": 8, "f": 8, "c": 8, "g": 8}


def _varint(x: int, out: bytearray) -> None:
    while x >= 128:
        out.append((x & 127) | 128)
        x >>= 7
    out.append(x)


def _write_gz(path: Path, data: bytes) -> int:
    path.write_bytes(gzip.compress(data, 9, mtime=0))
    return path.stat().st_size


def _block(gc: GeneCloud, gene_index: dict[str, int], prefixes: tuple[str, ...]):
    labels, codes, keys, buf = [], [], [], bytearray()
    for c in sorted(gc.members):
        if not c.startswith(prefixes):
            continue
        mem = gc.members[c]
        if len(mem) < 2:
            continue
        idx = sorted(gene_index[g] for g in mem)
        labels.append(gc.label[c])
        codes.append(LAYER_CODE.get(gc.kind[c], "g"))
        keys.append(c)
        _varint(len(idx), buf)
        prev = 0
        for i in idx:
            _varint(i - prev, buf)
            prev = i
    return labels, codes, keys, buf


def export(annotation: str | Path, outdir: str | Path, backgrounds: dict[str, list[str]] | None = None,
           phrase_min_genes: int = 5, phrase_min_npmi: float = 0.35) -> Path:
    """Write the web data files. ``backgrounds``: optional named gene lists flagged per gene (e.g. ATH1)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    engines = {}
    for preset in TEXT_PRESETS:
        engines[preset] = GeneCloud(annotation, layers=("word", "phrase", "keyword", "go"), text_fields=preset,
                                    phrase_min_genes=phrase_min_genes, phrase_min_npmi=phrase_min_npmi)
    ref = engines["all"]
    universe = sorted(set().union(*(e.concepts.keys() for e in engines.values())))
    gi = {g: i for i, g in enumerate(universe)}
    files, blocks = {}, {}

    def emit(name, eng, prefixes):
        labels, codes, keys, buf = _block(eng, gi, prefixes)
        size = _write_gz(outdir / f"index_{name}.bin.gz", bytes(buf))
        doc = {"labels": labels, "layers": "".join(codes)}
        if name == "go":                                   # concept keys are GO ids (labels are GO names)
            doc["ids"] = [k[2:] for k in keys]
        size += _write_gz(outdir / f"concepts_{name}.json.gz",
                          json.dumps(doc, ensure_ascii=False, separators=(",", ":")).encode())
        blocks[name] = {"concepts": len(labels), "bytes": size}

    for preset, eng in engines.items():
        emit(f"wp_{preset}", eng, ("w:", "p:"))
    emit("kw", ref, ("k:",))
    emit("go", ref, ("g:",))

    # per-gene metadata
    bgs = {k: set(ref.normalise(v)) for k, v in (backgrounds or {}).items()}
    sym, desc, flags, bits = [], [], [], {p: [] for p in TEXT_PRESETS}
    for g in universe:
        a = ref.genes.get(g)
        sym.append(a.symbol if a and a.symbol else "")
        d = ""
        if a:
            d = (a.tair_description[0] if a.tair_description else a.protein_names[0] if a.protein_names
                 else a.description) or ""
        desc.append(d[:140])
        f = 1 if (a and a.biotype == "protein_coding") else 0
        for j, name in enumerate(bgs):
            if g in bgs[name]:
                f |= 2 << j
        flags.append(f)
        for preset, eng in engines.items():
            b = 0
            for c in eng.concepts.get(g, ()):
                b |= LAYER_BIT.get(c[0], 0) if c[0] in "wpk" else 8
            bits[preset].append(b)
    _write_gz(outdir / "genes.json.gz", json.dumps(
        {"ids": universe, "symbols": sym, "desc": desc, "flags": flags, "bits": bits},
        ensure_ascii=False, separators=(",", ":")).encode())
    # full annotation texts, loaded by the web page only when a word is clicked (as in GeneCloud 2015, Fig. 1E)
    def cap(x, n):
        x = " ".join(x) if isinstance(x, list) else (x or "")
        return x if len(x) <= n else x[:n].rsplit(" ", 1)[0] + " …"
    texts = []
    for g in universe:
        a = ref.genes.get(g)
        if a is None:
            texts.append(["", "", "", "", "", ""])
            continue
        texts.append([cap(" | ".join(a.tair_description), 300), cap(a.tair_summary, 1500),
                      cap(" | ".join(a.tair_computational), 300),
                      cap(" | ".join(a.protein_names[:4]) + (" — " + a.function if a.function else ""), 1200),
                      cap(" | ".join(a.phenotypes), 800), a.biotype])
    _write_gz(outdir / "texts.json.gz", json.dumps(texts, ensure_ascii=False, separators=(",", ":")).encode())
    meta = {
        "genecloud_web": 1,
        "annotation": ref.meta,
        "n_genes": len(universe),
        "text_presets": {k: list(v) for k, v in TEXT_PRESETS.items()},
        "backgrounds": {"protein_coding": 1, **{name: 2 << j for j, name in enumerate(bgs)}},
        "background_sizes": {"protein_coding": sum(1 for x in flags if x & 1),
                             **{name: sum(1 for x in flags if x & (2 << j)) for j, name in enumerate(bgs)}},
        "blocks": blocks,
        "id_regex": ref.species.id_regex,
        "phrase_min_genes": phrase_min_genes, "phrase_min_npmi": phrase_min_npmi,
    }
    (outdir / "meta.json").write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    return outdir
