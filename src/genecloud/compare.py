"""Compare several gene sets (e.g. clusters): one dot plot of concepts x sets."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .cloud import CMAP_LIGHT


def compare_table(results: dict, per_set: int = 8) -> pd.DataFrame:
    """Long table (set, concept, fdr, fold, k) for the union of the top representatives of every set."""
    keep, disp = [], {}
    for r in results.values():
        top = r.significant.head(per_set)
        keep += list(top.label)
        for a, b in zip(top.label, top.display):
            disp.setdefault(a, b)
    keep = list(dict.fromkeys(keep))
    rows = []
    for name, r in results.items():
        t = r.table.drop_duplicates("label").set_index("label")
        for c in keep:
            if c in t.index:
                x = t.loc[c]
                rows.append(dict(set=name, concept=disp[c], label=c, fdr=x.fdr, p=x.p, fold=x.fold, k=x.k,
                                 significant=x.status == "significant"))
            else:
                rows.append(dict(set=name, concept=disp[c], label=c, fdr=1.0, p=1.0, fold=np.nan, k=0, significant=False))
    return pd.DataFrame(rows)


def draw_compare(results: dict, path, per_set: int = 8, title: str = "GeneCloud — comparison"):
    d = compare_table(results, per_set)
    if d.empty:
        raise ValueError("no significant concept in any set")
    concepts = list(dict.fromkeys(d.concept))
    sets = list(results)
    fig, ax = plt.subplots(figsize=(1.6 + 1.25 * len(sets), 1.4 + 0.27 * len(concepts)))
    lf = np.log2(d.fold.fillna(1).clip(lower=1))
    vmax = max(2, np.ceil(lf.max()))
    for _, r in d.iterrows():
        if r.k == 0:
            continue
        x, y = sets.index(r.set), concepts.index(r.concept)
        s = 25 + 60 * min(-np.log10(max(r.fdr, 1e-12)), 8)
        c = CMAP_LIGHT(min(np.log2(max(r.fold, 1)) / vmax, 1))
        ax.scatter(x, y, s=s, color=c if r.significant else "none", edgecolor=c, lw=1.4, zorder=3)
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels(sets, rotation=30, ha="right")
    ax.set_yticks(range(len(concepts)))
    ax.set_yticklabels(concepts)
    ax.invert_yaxis()
    ax.set_xlim(-0.6, len(sets) - 0.4)
    ax.grid(color="#eceff3", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, loc="left", fontweight="bold")
    for v in (0.05, 1e-3, 1e-6):
        ax.scatter([], [], s=25 + 60 * -np.log10(v), color="#9ca3af", label=f"FDR {v:g}")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, title="filled = FDR ≤ cut-off\nsize = −log10 FDR\ncolour = fold", title_fontsize=7.5)
    plt.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200 if str(path).endswith(".png") else None)
    plt.close(fig)
    return d
