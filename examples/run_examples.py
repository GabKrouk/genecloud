"""Re-run the examples of the GeneCloud paper (Krouk et al., 2015, Mol. Plant 8:971) with GeneCloud 2.

    genecloud build --species arabidopsis --out arabidopsis.jsonl.gz      # once
    python examples/run_examples.py arabidopsis.jsonl.gz                  # -> examples/results/

All lists come from ATH1 microarray experiments, so the background is every gene on the ATH1 array
(examples/lists/background_ATH1.txt). The lists are rebuilt from the public sources by prepare_lists.py.
"""
import sys
from pathlib import Path

from genecloud import GeneCloud
from genecloud.compare import draw_compare

HERE = Path(__file__).parent
LISTS, OUT, IMG = HERE / "lists", HERE / "results", HERE / "images"
HORMONES = ["IAA", "BL", "CK", "GA", "ABA", "ACC", "MeJA"]

EXAMPLES = {
    "hormone_IAA_up": "Auxin (IAA)-induced genes · Nemhauser et al. 2006",
    "PHR1_induced": "Genes induced by PHR1 over-expression · Bustos et al. 2010",
    "nitrate_responsive": "Nitrate-responsive genes · Canales et al. 2014",
}


def ids(path):
    return [l.strip() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]


def main(annotation):
    OUT.mkdir(exist_ok=True)
    IMG.mkdir(exist_ok=True)
    gc = GeneCloud(annotation, layers=("word", "phrase", "keyword"))
    bg = ids(LISTS / "background_ATH1.txt")

    summary = []
    for name, title in EXAMPLES.items():
        res = gc.run(ids(LISTS / f"{name}.txt"), background=bg)
        res.to_tsv(OUT / f"{name}.tsv")
        gc.html(res, OUT / f"{name}.html", title=title)
        gc.plot(res, OUT / f"{name}.pdf", title=title)
        gc.plot(res, IMG / f"{name}.png", title=title)
        summary.append((name, res))

    # the seven hormones of Nemhauser et al. 2006, induced genes, side by side
    hormones = {}
    for h in HORMONES:
        res = gc.run(ids(LISTS / f"hormone_{h}_up.txt"), background=bg)
        res.to_tsv(OUT / f"hormone_{h}_up.tsv")
        hormones[f"{h} up"] = res
    d = draw_compare(hormones, IMG / "hormones_compare.png", per_set=6,
                     title="Genes induced by seven hormones (Nemhauser et al. 2006)")
    d.to_csv(OUT / "hormones_compare.tsv", sep="\t", index=False)

    for name, res in summary + list(hormones.items()):
        sig = res.significant
        top = ", ".join(sig.display.head(8)) if len(sig) else "-"
        print(f"{name:22s} n={len(res.study):4d}/{len(res.study) + len(res.missing):4d}  "
              f"significant={len(sig):3d}  top: {top}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "arabidopsis.jsonl.gz")
