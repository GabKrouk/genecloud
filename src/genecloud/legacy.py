"""The original GeneCloud test (Krouk, 2013), kept for reproducibility and comparison.

Word frequencies (all occurrences, not genes) in the set's descriptions are compared with ``sampling`` random
gene lists of the same size drawn from the genome; a word is called significant if no random list reaches
its frequency. The enrichment ratio is (freq / n) / (genome count / genome size), as in the R package.
"""
from __future__ import annotations

import re
from collections import Counter

import numpy as np
import pandas as pd

_PUNCT = re.compile(r"[!\"#$%&'()*+,\-./:;<=>?@\[\\\]^_`{|}~]")
_MEANINGLESS = {"protein", "proteincoding", "expressed", "proteins", "plants", "blast", "has", "source", "blink",
                "hits", "ncbi", "match", "contains", "best", "interpro", "involved", "functions", "plant", "encodes"}


def _freq(texts, stop):
    c = Counter()
    for t in texts:
        t = _PUNCT.sub("", t)
        c.update(w for w in t.lower().split() if len(w) >= 3 and w not in stop)
    return c


def legacy_genecloud(gc, genes, sampling: int = 100, seed: int = 1) -> pd.DataFrame:
    from .text import ENGLISH_STOPWORDS
    rng = np.random.default_rng(seed)
    text = {g: " ".join(gc._text(a, ("names", "description", "protein_names", "function", "families")))
            for g, a in gc.genes.items()}
    universe = np.array(sorted(g for g, t in text.items() if t))
    study = [g for g in gc.normalise(genes) if g in set(universe)]
    obs = _freq((text[g] for g in study), ENGLISH_STOPWORDS)
    genome = _freq(text.values(), ENGLISH_STOPWORDS)
    rdm = Counter()
    for _ in range(sampling):
        r = _freq((text[g] for g in rng.choice(universe, len(study), replace=False)), ENGLISH_STOPWORDS)
        for w, m in obs.items():
            if r.get(w, 0) >= m:
                rdm[w] += 1
    t = pd.DataFrame({"word": list(obs), "freq": [obs[w] for w in obs], "rdm": [rdm[w] for w in obs],
                      "genome_count": [genome[w] for w in obs]})
    t["enrich_ratio"] = (t.freq / len(study)) / (t.genome_count / len(universe))
    t["significant"] = (t.rdm == 0) & ~t.word.isin(_MEANINGLESS)
    return t.sort_values(["significant", "freq"], ascending=False).reset_index(drop=True)
