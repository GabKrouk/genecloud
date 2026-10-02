"""GeneCloud engine: concept index over the genome, enrichment test, redundancy reduction."""
from __future__ import annotations

import gzip
import hashlib
import math
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import hypergeom

from . import annotation as ann
from .species import SPECIES
from .text import Tokenizer, detect_phrases, phrases_in

LAYERS = ("word", "phrase", "keyword", "go")
# which annotation texts feed the word / phrase layers
TEXT_PRESETS = {
    # everything: TAIR (descriptions, curator summaries, phenotypes) + NCBI + UniProt + GO term names
    "all": ("names", "description", "tair_description", "tair_summary", "tair_computational", "phenotypes",
            "protein_names", "function", "families", "go_names"),
    # TAIR only, the closest to GeneCloud 2013/2015 (which used TAIR descriptions), but with today's text
    "tair": ("names", "tair_description", "tair_summary", "tair_computational", "phenotypes"),
    # the GeneCloud 2.0 set (no TAIR)
    "open": ("names", "description", "protein_names", "function", "families", "go_names"),
}
TEXT_FIELDS = TEXT_PRESETS["all"]
GO_NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}


def _stem(w: str) -> str:
    for suf in ("ers", "er", "ions", "ion", "ing", "ed"):
        if w.endswith(suf) and len(w) - len(suf) >= 4:
            return w[: -len(suf)]
    return w


def log10_hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """log10 P(X >= k), X ~ hypergeometric(N, K, n); exact even far below 1e-300."""
    p = float(hypergeom.sf(k - 1, N, K, n))
    if p > 1e-280:
        return math.log10(p) if p > 0 else -math.inf
    i = np.arange(k, min(K, n) + 1)
    return float(logsumexp(hypergeom.logpmf(i, N, K, n)) / math.log(10))


def log10_adjust(log10_p: np.ndarray, m: int | None = None, method: str = "BH") -> np.ndarray:
    """Benjamini-Hochberg (or -Yekutieli) adjusted P values, in log10, for ``m`` hypotheses (>= len(p))."""
    lp = np.asarray(log10_p, float)
    n = len(lp)
    if n == 0:
        return lp
    m = n if m is None else max(m, n)
    c = math.log10(sum(1.0 / i for i in range(1, m + 1))) if method == "BY" else 0.0
    o = np.argsort(lp, kind="stable")
    q = lp[o] + math.log10(m) - np.log10(np.arange(1, n + 1)) + c
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 0.0)
    return out


def bh(p: np.ndarray, m: int | None = None) -> np.ndarray:
    """Benjamini-Hochberg adjusted P values. ``m`` = total number of hypotheses (>= len(p)); hypotheses
    that are not listed are taken as P = 1 (they rank last and do not change the adjusted values above)."""
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    m = n if m is None else max(m, n)
    o = np.argsort(p, kind="stable")
    q = p[o] * m / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


@dataclass
class CloudResult:
    table: pd.DataFrame                  # every tested concept
    study: list                          # study genes used (annotated, in background)
    missing: list                        # study genes not found / not in background
    background_size: int
    params: dict = field(default_factory=dict)
    not_in_background: list = field(default_factory=list)
    not_annotated: list = field(default_factory=list)      # in the background but without annotation

    @property
    def significant(self) -> pd.DataFrame:
        t = self.table
        return t[(t.status == "significant") & t.representative]

    @property
    def trend(self) -> pd.DataFrame:
        t = self.table
        return t[(t.status == "trend") & t.representative]

    def to_tsv(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.table.to_csv(path, sep="\t", index=False)


def _code_tag() -> str:
    """Changes whenever the text-processing code changes, so stale concept caches are rebuilt."""
    src = b"".join((Path(__file__).parent / f).read_bytes() for f in ("text.py", "core.py"))
    return hashlib.md5(src).hexdigest()[:8]


class GeneCloud:
    """Build once per annotation, then call :meth:`run` on any gene list.

    >>> gc = GeneCloud("arabidopsis.jsonl.gz")
    >>> res = gc.run(my_genes, background=expressed_genes)
    >>> gc.plot(res, "cloud.pdf")
    """

    def __init__(self, annotation_path: str | Path, species: str = "arabidopsis",
                 layers: Sequence[str] = ("word", "phrase", "keyword"),
                 text_fields: Sequence[str] | str = TEXT_FIELDS, keep_symbols: bool = False,
                 phrase_min_genes: int = 5, phrase_min_npmi: float = 0.35, cache: bool = True):
        self.annotation_path = Path(annotation_path)
        if isinstance(text_fields, str):
            text_fields = TEXT_PRESETS[text_fields]
        self.text_fields = tuple(text_fields)
        self.species = SPECIES[species] if isinstance(species, str) else species
        self.layers = tuple(layers)
        bad = set(self.layers) - set(LAYERS)
        if bad:
            raise ValueError(f"unknown layer(s) {bad}; choose from {LAYERS}")
        self.genes, self.meta = ann.load(self.annotation_path)
        self.onto = ann.ontology_for(self.annotation_path)
        key = repr((self.layers, tuple(text_fields), keep_symbols, phrase_min_genes, phrase_min_npmi,
                    self.annotation_path.stat().st_size, self.meta.get("built"), _code_tag()))
        cache_path = self.annotation_path.with_suffix(".idx." + hashlib.md5(key.encode()).hexdigest()[:10] + ".pkl.gz")
        if cache and cache_path.exists():
            with gzip.open(cache_path, "rb") as fh:
                self.concepts, self.kind, self.label = pickle.load(fh)
        else:
            self._index(text_fields, keep_symbols, phrase_min_genes, phrase_min_npmi)
            if cache:
                with gzip.open(cache_path, "wb") as fh:
                    pickle.dump((self.concepts, self.kind, self.label), fh)
        # inverted index
        self.members: dict[str, set] = {}
        for g, cs in self.concepts.items():
            for c in cs:
                self.members.setdefault(c, set()).add(g)

    # ------------------------------------------------------------------ index
    def _text(self, g: ann.Gene, fields: Sequence[str]) -> list[str]:
        out = []
        for f in fields:
            if f == "go_names":
                if self.onto is not None:
                    out += [self.onto.name(t) for t in g.go]
            else:
                v = getattr(g, f)
                out += v if isinstance(v, list) else [v]
        return [t for t in out if t]

    def _index(self, fields, keep_symbols, ph_min, ph_npmi):
        tok = Tokenizer(self.species.id_pattern, self.species.extra_stopwords, keep_symbols)
        runs = {gid: [r for t in self._text(g, fields) for r in tok.sentences(t)] for gid, g in self.genes.items()}
        phrases = detect_phrases(runs, ph_min, ph_npmi) if "phrase" in self.layers else set()
        self.concepts, self.kind, self.label = {}, {}, {}
        for gid, g in self.genes.items():
            cs = set()
            if "word" in self.layers:
                for w in {w for r in runs[gid] for w in r}:
                    cs.add("w:" + w)
                    self.kind["w:" + w], self.label["w:" + w] = "word", w
            if "phrase" in self.layers:
                for p in phrases_in(runs[gid], phrases):
                    cs.add("p:" + p)
                    self.kind["p:" + p], self.label["p:" + p] = "phrase", p
            if "keyword" in self.layers:
                for k in g.keywords:
                    if k in ("Reference proteome", "Proteomics identification", "3D-structure", "Direct protein sequencing",
                             "Alternative splicing", "Signal", "Phosphoprotein", "Acetylation"):
                        continue
                    cs.add("k:" + k)
                    self.kind["k:" + k], self.label["k:" + k] = "keyword", k
            if "go" in self.layers and self.onto is not None:
                for t in {a for d in g.go for a in self.onto.ancestors(d)}:
                    if t in ("GO:0008150", "GO:0003674", "GO:0005575"):
                        continue
                    cs.add("g:" + t)
                    self.kind["g:" + t] = "GO " + GO_NS.get(self.onto.namespace(t), "")
                    self.label["g:" + t] = self.onto.name(t)
            if cs:
                self.concepts[gid] = cs

    # ------------------------------------------------------------------ helpers
    def normalise(self, ids: Iterable[str]) -> list[str]:
        out = []
        for x in ids:
            x = str(x).strip()
            m = self.species.id_pattern.search(x)
            if m:
                out.append(self.species.normalise(m.group(0)))
        return list(dict.fromkeys(out))

    def symbol(self, gid: str) -> str:
        g = self.genes.get(gid)
        return (g.symbol if g and g.symbol else gid)

    # ------------------------------------------------------------------ test
    def run(self, genes: Iterable[str], background: Iterable[str], min_genes: int = 2,
            fdr: float = 0.05, max_concept_frac: float = 0.25, merge_jaccard: float = 0.75,
            min_fold: float = 1.0, trend_p: float = 0.01, adjust: str = "BH") -> CloudResult:
        """Hypergeometric over-representation of every concept in ``genes`` vs ``background``.

        background       REQUIRED. The genes that could have been in the list (e.g. all genes detected in the
                         experiment, or all genes on the array). Study genes absent from it are ignored.

        min_genes        concept must be present in at least this many study genes (kills gene-name noise)
        max_concept_frac concepts carried by more than this fraction of the background are not tested
        merge_jaccard    significant concepts carried by (almost) the same study genes are merged; the most
                         significant one represents the group (e.g. 'ammonium' + 'ammonium transport')
        adjust           "BH" (Benjamini-Hochberg, default) or "BY" (Benjamini-Yekutieli, valid under any
                         dependence between concepts, more conservative)
        """
        if adjust not in ("BH", "BY"):
            raise ValueError("adjust must be 'BH' or 'BY'")
        if background is None or isinstance(background, str):
            raise ValueError("a background gene list is required (the genes that could have been in the list, "
                             "e.g. all expressed genes or all genes on the array)")
        req = self.normalise(genes)
        bg = set(self.normalise(background))
        if not bg:
            raise ValueError("the background contains no recognised gene identifier")
        universe = set(self.concepts) & bg
        study = [g for g in req if g in universe]
        missing = [g for g in req if g not in universe]
        not_in_background = [g for g in req if g not in bg]
        not_annotated = [g for g in missing if g in bg]
        N, n = len(universe), len(study)
        if n == 0:
            raise ValueError("none of the genes are annotated / in the background")
        counts_bg: dict[str, int] = {}
        for c, mem in self.members.items():
            k = len(mem & universe)
            if k:
                counts_bg[c] = k
        # Hypotheses = every concept that could reach min_genes in a list drawn from this background and is not
        # over-general: min_genes <= K <= max_concept_frac * N. This filter depends on the background only
        # (never on the study list), so Benjamini-Hochberg over this family controls the FDR; concepts absent
        # from the study list are part of the family (P = 1).
        testable = {c for c, K in counts_bg.items() if min_genes <= K <= max_concept_frac * N}
        hits: dict[str, list] = {}
        for g in study:
            for c in self.concepts[g]:
                if c in testable:
                    hits.setdefault(c, []).append(g)
        rows = []
        for c, gs in hits.items():
            K, k = counts_bg[c], len(gs)
            lp = log10_hypergeom_sf(k, N, K, n)
            rows.append((c, self.kind[c], self.label[c], k, n, K, N, (k / n) / (K / N), lp, sorted(gs)))
        t = pd.DataFrame(rows, columns=["concept", "layer", "label", "k", "n", "K", "N", "fold", "log10_p", "genes"])
        m_tests = len(testable)
        t["tested"] = t.k >= min_genes
        # concepts below min_genes can never be called: they enter the correction with P = 1 (textbook-valid)
        lp_adj = np.where(t.tested, t.log10_p, 0.0) if len(t) else np.array([])
        t["log10_fdr"] = log10_adjust(lp_adj, m_tests, adjust) if len(t) else []
        t["p"] = 10.0 ** t.log10_p
        t["fdr"] = 10.0 ** t.log10_fdr
        t = t.sort_values(["log10_fdr", "log10_p", "k", "concept"], ascending=[True, True, False, True],
                          kind="mergesort").reset_index(drop=True)
        # redundancy: greedy grouping of candidate concepts (significant or trend) in order of significance.
        # i joins an earlier representative j if they are carried by (almost) the same study genes, or if the
        # words of one are contained in the other and they share at least half of the genes of the larger one
        # ('high' -> 'high affinity'; but 'phosphate' is not absorbed by the smaller 'pentose phosphate').
        t["status"] = np.where(t.tested & (t.fdr <= fdr) & (t.fold >= min_fold), "significant",
                               np.where(t.tested & (t.p <= trend_p) & (t.fold >= min_fold), "trend", ""))
        t["representative"] = False
        t["group"] = ""
        t["display"] = t.label
        cand = t.index[t.status != ""]
        reps: list[tuple[int, set, set]] = []
        for i in cand:
            gs, toks = set(t.at[i, "genes"]), {_stem(x) for x in str(t.at[i, "label"]).lower().split()}
            for j, rs, rt in reps:
                inter = len(gs & rs)
                jac = inter / len(gs | rs)
                contained = toks <= rt or rt <= toks
                if (jac >= merge_jaccard or (contained and inter / max(len(gs), len(rs)) >= 0.5)) \
                        and (t.at[i, "status"] == t.at[j, "status"] or t.at[j, "status"] == "significant"):
                    t.at[i, "group"] = t.at[j, "concept"]
                    break
            else:
                reps.append((i, gs, toks))
                t.at[i, "representative"] = True
                t.at[i, "group"] = t.at[i, "concept"]
        # a lone word is shown with the phrase that (nearly) all its genes share ('distance' -> 'long distance'):
        # >= 90 % of the genes, so that the statistics drawn for the word are those of what is written
        for i in t.index[t.representative & (t.layer == "word")]:
            w = t.at[i, "label"]
            gs = t.at[i, "genes"]
            cnt: dict[str, int] = {}
            for g in gs:
                for c in self.concepts[g]:
                    if c.startswith("p:") and w in c[2:].split():
                        cnt[c] = cnt.get(c, 0) + 1
            good = [c for c, v in cnt.items() if v >= max(2, math.ceil(0.9 * len(gs)))]
            if good:   # most genes, then the most specific phrase, then alphabetical (deterministic)
                best = min(good, key=lambda c: (-cnt[c], len(self.members[c]), c))
                t.at[i, "display"] = self.label[best]
        members = {}
        for i in cand:
            members.setdefault(t.at[i, "group"], []).append(t.at[i, "label"])
        t["merged"] = [", ".join(members.get(t.at[i, "concept"], [])[1:]) if t.at[i, "representative"] else ""
                       for i in t.index]
        glabel = dict(zip(t.concept, t.label))
        t["group"] = [glabel.get(g, "") for g in t.group]
        for c in ("label", "display", "group", "merged"):
            t[c] = t[c].astype(str).str.replace("\u2011", "-")
        t["symbols"] = [", ".join(self.symbol(g) for g in gs) for gs in t.genes]
        t["genes"] = [", ".join(gs) for gs in t.genes]
        params = dict(min_genes=min_genes, fdr=fdr, max_concept_frac=max_concept_frac,
                      merge_jaccard=merge_jaccard, min_fold=min_fold, trend_p=trend_p, adjust=adjust, layers=self.layers,
                      text_fields=self.text_fields, n_hypotheses=m_tests,
                      background_genes=len(bg), background_annotated=N,
                      annotation=self.meta)
        cols = ["concept", "layer", "label", "k", "n", "K", "N", "fold", "p", "fdr", "log10_p", "log10_fdr", "tested",
                "status", "representative", "group", "display", "merged", "genes", "symbols"]
        return CloudResult(t[cols], study, missing, N, params, not_in_background, not_annotated)

    # ------------------------------------------------------------------ convenience
    def plot(self, result: CloudResult, path: str | Path, **kw):
        from .cloud import draw_cloud
        return draw_cloud(result, path, **kw)

    def html(self, result: CloudResult, path: str | Path, **kw):
        from .report import write_html
        return write_html(result, path, gc=self, **kw)
