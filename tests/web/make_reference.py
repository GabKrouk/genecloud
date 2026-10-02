"""Reference results from the Python engine for tests/web/test_engine.mjs."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

from genecloud import GeneCloud
from genecloud.core import log10_hypergeom_sf

ann, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
L = Path("examples/lists")
ids = lambda f: [l.strip() for l in open(f) if l.strip() and not l.startswith("#")]
CASES = [
    ("IAA_all", "hormone_IAA_up", "background_ATH1", dict(preset="all", layers=["word", "phrase", "keyword"])),
    ("PHR1_all", "PHR1_induced", "background_ATH1", dict(preset="all", layers=["word", "phrase", "keyword"])),
    ("nitrate_all", "nitrate_responsive", "background_ATH1", dict(preset="all", layers=["word", "phrase", "keyword"])),
    ("GA_tair_go", "hormone_GA_up", "background_ATH1", dict(preset="tair", layers=["word", "phrase", "go"])),
    ("CK_open_strict", "hormone_CK_up", "background_ATH1", dict(preset="open", layers=["word", "keyword", "go"],
                                                                 minGenes=3, fdr=0.01, maxFrac=0.1)),
    ("MeJA_all_BY", "hormone_MeJA_up", "background_ATH1", dict(preset="all", layers=["word", "phrase", "keyword", "go"],
                                                                adjust="BY")),
]
cases = []
for name, genes, bg, o in CASES:
    gc = GeneCloud(ann, layers=tuple(o["layers"]), text_fields=o["preset"])
    kw = dict(min_genes=o.get("minGenes", 2), fdr=o.get("fdr", 0.05), max_concept_frac=o.get("maxFrac", 0.25),
              adjust=o.get("adjust", "BH"))
    res = gc.run(ids(L / f"{genes}.txt"), ids(L / f"{bg}.txt"), **kw)
    res.table.to_csv(out / f"{name}.tsv", sep="\t", index=False)
    cases.append(dict(name=name, genes=str(L / f"{genes}.txt"), background=str(L / f"{bg}.txt"), opts=o))
    print(name, len(res.table), len(res.significant))
(out / "cases.json").write_text(json.dumps(cases, indent=1))
rng = np.random.default_rng(1)
hg = []
for _ in range(3000):
    N = int(rng.integers(20, 40000)); K = int(rng.integers(1, N)); n = int(rng.integers(1, N))
    lo, hi = max(0, n - (N - K)), min(K, n)
    k = int(rng.integers(lo, hi + 1))
    hg.append([k, N, K, n, float(hypergeom.sf(k - 1, N, K, n)), log10_hypergeom_sf(k, N, K, n)])
for k, N, K, n in [(400, 30000, 500, 600), (300, 21633, 1200, 2136), (90, 21633, 193, 2136), (60, 40000, 60, 60)]:
    hg.append([k, N, K, n, float(hypergeom.sf(k - 1, N, K, n)), log10_hypergeom_sf(k, N, K, n)])
(out / "hypergeom.json").write_text(json.dumps(hg))
