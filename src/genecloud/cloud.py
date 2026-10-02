"""Word-cloud rendering: word size = statistical enrichment (-log10 FDR), colour = fold enrichment.

The layout is deterministic (no random placement): concepts are placed from the most to the least
significant along an Archimedean spiral with exact text bounding boxes, so two runs give the same picture.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex

# perceptually ordered, colour-blind friendly ramps
CMAP_LIGHT = LinearSegmentedColormap.from_list("gc_light", ["#8fc9c4", "#2f8f9d", "#2a5ea8", "#5b3a9e", "#b0216f"])
CMAP_DARK = LinearSegmentedColormap.from_list("gc_dark", ["#6b7280", "#38bdf8", "#facc15", "#f97316", "#ef4444"])
LAYER_STYLE = {"keyword": "italic"}
THEMES = {"light": dict(bg="#ffffff", fg="#111827", sub="#6b7280", trend="#c7cbd1", cmap=CMAP_LIGHT),
          "dark": dict(bg="#0b0f14", fg="#f3f4f6", sub="#9ca3af", trend="#4b5563", cmap=CMAP_DARK)}


@dataclass
class Placed:
    text: str
    x: float          # centre, points
    y: float
    w: float
    h: float
    size: float
    color: str
    style: str
    weight: str
    row: dict         # the result row (concept statistics)


def _spiral(items, width, height, fig, renderer, pad=2.0):
    placed, out = [], []
    cx, cy = width / 2, height / 2
    ax_ratio = (width / height) ** 0.5
    for text, fs, kw in items:
        t = fig.text(0, 0, text, fontsize=fs, **kw)
        bb = t.get_window_extent(renderer)
        t.remove()
        w, h = bb.width * 72 / fig.dpi + pad, bb.height * 72 / fig.dpi + pad
        theta, pos = 0.0, None
        while theta < 400:
            r = 2.2 * theta
            x = cx + r * math.cos(theta) * ax_ratio - w / 2
            y = cy + r * math.sin(theta) / ax_ratio - h / 2
            if 0 <= x and x + w <= width and 0 <= y and y + h <= height and \
                    all(x + w <= a or a + c <= x or y + h <= b or b + d <= y for a, b, c, d in placed):
                pos = (x, y, w, h)
                break
            theta += 0.35 / (1 + theta / 60)
        if pos:
            placed.append(pos)
        out.append(pos)
    return out


def layout(result, width=792, height=450, top=60, show_trend=True, max_trend=25, theme="light",
           max_font=44, min_font=9, layers=None, size_by="fdr"):
    """Compute word positions (in points). Returns (placed words, fold-enrichment Normalize, theme dict)."""
    th = THEMES[theme]
    sig = result.significant
    tr = result.trend if show_trend else result.trend.iloc[:0]
    if layers:
        sig, tr = sig[sig.layer.isin(layers)], tr[tr.layer.isin(layers)]
    sig, tr = sig.head(top), tr.head(max_trend)
    lf = np.log2(sig.fold.values) if len(sig) else np.array([1.0])
    lo = max(0.0, np.floor(lf.min()))
    norm = Normalize(vmin=lo, vmax=max(lo + 1.0, np.ceil(lf.max())))
    norm.size_key = []
    if not len(sig) and not len(tr):
        return [], norm, th
    fig = plt.figure(figsize=(width / 72, height / 72))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    by_fold = size_by == "fold"
    s = np.log2(sig.fold.values) if by_fold else -sig.log10_fdr.values
    scale = 1.0
    for _ in range(10):                                   # shrink until every significant concept fits
        hi = max(min_font * 1.2, max_font * scale)
        if len(s):
            smin = s.min() if by_fold else min(s.min(), -math.log10(result.params["fdr"]))
            smax = s.max()
            fs = min_font + (hi - min_font) * ((s - smin) / (smax - smin) if smax > smin else np.ones_like(s))
        else:
            fs = np.array([])
        items, meta = [], []
        for (_, r), f, l in zip(sig.iterrows(), fs, lf if len(sig) else []):
            kw = dict(fontweight="bold", style=LAYER_STYLE.get(r.layer, "normal"), family="DejaVu Sans")
            items.append((r.display, float(f), kw))
            meta.append((to_hex(th["cmap"](norm(l))), r))
        for _, r in tr.iterrows():
            kw = dict(fontweight="normal", style=LAYER_STYLE.get(r.layer, "normal"), family="DejaVu Sans")
            items.append((r.display, min_font * 0.95, kw))
            meta.append((th["trend"], r))
        pos = _spiral(items, width, height, fig, rend)
        if all(p is not None for p in pos[:len(sig)]):
            break
        scale *= 0.85
    plt.close(fig)
    if len(s):                                            # size key: a few round FDR values and their font size
        def size_of(v):
            return min_font + (hi - min_font) * ((v - smin) / (smax - smin) if smax > smin else 1.0)
        if by_fold:
            marks = list(range(int(math.ceil(smin)), int(math.floor(smax)) + 1)) or [smin]
            if len(marks) > 4:
                step = math.ceil(len(marks) / 4)
                marks = [v for i, v in enumerate(marks) if i % step == 0 or i == len(marks) - 1]
            norm.size_key = [(2.0 ** v, size_of(v), f"×{2 ** v:.2g}") for v in marks]
        else:
            cut = -math.log10(result.params["fdr"])
            top_mark = int(math.floor(smax))
            marks = [smin] + ([round((smin + top_mark) / 2)] if top_mark - smin > 4 else []) + \
                    ([top_mark] if top_mark - smin > 1.5 else [])
            norm.size_key = [(10.0 ** -v, size_of(v), (f"{result.params['fdr']:g}" if abs(v - cut) < 1e-9
                                                       else f"10⁻{_sup(int(v))}")) for v in marks]
        norm.size_by = size_by
    out = []
    for (text, f, kw), p, (col, r) in zip(items, pos, meta):
        if p is not None:
            x, y, w, h = p
            out.append(Placed(text, x + w / 2, y + h / 2, w, h, f, col, kw["style"], kw["fontweight"], r.to_dict()))
    return out, norm, th


def _sup(n: int) -> str:
    return str(n).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))


def draw_cloud(result, path, top: int = 60, title: str | None = None, theme: str = "light",
               max_font: float = 44, min_font: float = 9, size: tuple = (11, 7.5), show_legend: bool = True,
               show_trend: bool = True, max_trend: int = 25, layers: tuple | None = None, size_by: str = "fdr"):
    """Draw the cloud of significant, non-redundant concepts (plus, in grey, nominal trends).

    Word size is proportional to -log10(FDR); colour encodes fold enrichment. ``path``: .pdf, .png or .svg.
    """
    W, H = size[0] * 72, (size[1] - 1.9) * 72       # cloud area between the header and the legends
    words, norm, th = layout(result, W, H, top, show_trend, max_trend, theme, max_font, min_font, layers, size_by)
    fig = plt.figure(figsize=size, facecolor=th["bg"])
    ax = fig.add_axes([0, 0.92 / size[1], 1, (size[1] - 1.9) / size[1]])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    for p in words:
        ax.text(p.x, p.y, p.text, fontsize=p.size, ha="center", va="center", color=p.color,
                fontweight=p.weight, style=p.style, family="DejaVu Sans")
    n_sig = sum(p.row["status"] == "significant" for p in words)
    if not n_sig:
        ax.text(W / 2, 0.06 * H, f"no concept passes FDR {result.params['fdr']}", ha="center", color=th["sub"], fontsize=11)
    n, N = len(result.study), result.background_size
    fig.text(0.02, 1 - 0.32 / size[1], title or "GeneCloud", fontsize=17, fontweight="bold", color=th["fg"], va="top")
    fig.text(0.02, 1 - 0.68 / size[1],
             f"{n} genes vs {N:,} background genes  ·  size = {'log2 fold enrichment' if size_by == 'fold' else '−log10 FDR (hypergeometric, BH)'}  ·  colour = fold enrichment  ·  "
             f"FDR ≤ {result.params['fdr']}, ≥ {result.params['min_genes']} genes, redundant concepts merged",
             fontsize=8.5, color=th["sub"], va="top")
    if show_legend and n_sig:
        cax = fig.add_axes([0.70, 0.40 / size[1], 0.26, 0.11 / size[1]])
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=th["cmap"]), cax=cax, orientation="horizontal")
        ticks = np.arange(np.ceil(norm.vmin), np.floor(norm.vmax) + 1)
        cb.set_ticks(ticks)
        cb.set_ticklabels([f"×{2 ** t:g}" for t in ticks])
        cb.ax.tick_params(labelsize=7.5, colors=th["sub"])
        cb.outline.set_edgecolor(th["sub"])
        cb.set_label("fold enrichment (colour)", fontsize=8, color=th["sub"])
        cb.ax.xaxis.set_label_position("top")
    if show_legend and n_sig and getattr(norm, "size_key", None):   # word size key (FDR)
        fig.text(0.02, 0.62 / size[1], "word size (log2 fold):" if size_by == "fold" else "word size (−log10 FDR):",
                 fontsize=7.5, color=th["sub"], va="bottom")
        xpt = 0.02 * size[0] * 72 + 100
        for v, f, lab in norm.size_key:
            t = fig.text(xpt / (size[0] * 72), 0.62 / size[1], lab, fontsize=min(f, 15), color=th["sub"],
                         fontweight="bold", va="bottom")
            xpt += len(lab) * min(f, 15) * 0.62 + 14
    has_trend = any(p.row["status"] == "trend" for p in words)
    fig.text(0.02, 0.36 / size[1], "italic = UniProt keyword" +
             (f";  grey = nominal trend (uncorrected P ≤ {result.params.get('trend_p')}, not significant after correction)"
              if has_trend else ""),
             fontsize=7.5, color=th["sub"])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, facecolor=th["bg"], dpi=200 if str(path).endswith(".png") else None)
    plt.close(fig)
    return path
