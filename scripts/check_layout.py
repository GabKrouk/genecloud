"""Visual QA of the clouds: no two texts overlap, nothing leaves the figure, nothing sits on the colour bar.

    python scripts/check_layout.py arabidopsis.jsonl.gz
"""
import itertools
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

from genecloud import GeneCloud
from genecloud import cloud as C

ids = lambda f: [l.strip() for l in open(f) if l.strip() and not l.startswith("#")]
gc = GeneCloud(sys.argv[1])
bg = ids("examples/lists/background_ATH1.txt")
orig, bad = plt.Figure.savefig, 0
for name in ["hormone_IAA_up", "PHR1_induced", "nitrate_responsive", "hormone_CK_up", "hormone_GA_up"]:
    res = gc.run(ids(f"examples/lists/{name}.txt"), bg)
    for size_by in ("fdr", "fold"):
        for theme in ("light", "dark"):
            cap = {}
            plt.Figure.savefig = lambda self, *a, **k: (cap.__setitem__("f", self), orig(self, *a, **k))[1]
            C.draw_cloud(res, "/tmp/genecloud_check.png", size_by=size_by, theme=theme, title=name)
            plt.Figure.savefig = orig
            fig = cap["f"]
            cv = FigureCanvasAgg(fig)
            cv.draw()
            rd = cv.get_renderer()
            texts = list(fig.axes[0].texts) + list(fig.texts)
            for a in fig.axes[1:]:
                texts += list(a.get_xticklabels()) + [a.xaxis.label]
            boxes = [t.get_window_extent(rd) for t in texts if t.get_text()]
            ov = sum(1 for a, b in itertools.combinations(boxes, 2) if a.overlaps(b))
            out = sum(1 for b in boxes if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.x1 + 1 or b.y1 > fig.bbox.y1 + 1)
            cb = [a.get_window_extent(rd) for a in fig.axes[1:]]
            on_cb = sum(1 for t in fig.axes[0].texts for c in cb if t.get_window_extent(rd).overlaps(c))
            shown = {t.get_text() for t in fig.axes[0].texts}
            lost = sum(1 for d in res.significant.display.head(60) if d not in shown)
            ok = not (ov or out or on_cb or lost)
            bad += not ok
            print(f"{'ok  ' if ok else 'FAIL'} {name:20s} {size_by:4s} {theme:5s} texts={len(boxes):3d} overlaps={ov} "
                  f"outside={out} on-colourbar={on_cb} significant-not-drawn={lost}")
sys.exit(1 if bad else 0)
