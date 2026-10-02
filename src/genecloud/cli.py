"""Command line interface:  genecloud {build,run,compare,legacy} ..."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__


def _read_ids(path: str) -> list[str]:
    text = sys.stdin.read() if path == "-" else Path(path).read_text()
    return [tok for line in text.splitlines() if not line.lstrip().startswith("#")
            for tok in line.replace(",", " ").replace("\t", " ").split()]


def _read_sets(path: str) -> dict[str, list[str]]:
    """Two-column file (gene, set), tab or comma separated, header optional."""
    sets: dict[str, list[str]] = {}
    with open(path, newline="") as fh:
        sample = fh.read(2048)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;") if sample else csv.excel_tab
        for row in csv.reader(fh, dialect):
            if len(row) >= 2 and row[0].strip():
                sets.setdefault(row[1].strip(), []).append(row[0].strip())
    return sets


def _engine(args):
    from .core import GeneCloud
    return GeneCloud(args.annotation, species=args.species, layers=tuple(args.layers.split(",")))


def _run_opts(args):
    return dict(min_genes=args.min_genes, fdr=args.fdr, max_concept_frac=args.max_frac,
                merge_jaccard=args.merge, trend_p=args.trend_p)


def cmd_build(args):
    from .annotation import build, download
    src = Path(args.sources)
    if not args.offline:
        download(args.species, src)
    out = build(src, args.species, args.out)
    print(out)


def cmd_run(args):
    gc = _engine(args)
    bg = _read_ids(args.background)
    res = gc.run(_read_ids(args.genes), background=bg, **_run_opts(args))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    title = args.title or f"GeneCloud · {Path(args.genes).stem}"
    res.to_tsv(f"{out}.tsv")
    for fmt in args.formats.split(","):
        if fmt == "html":
            gc.html(res, f"{out}.html", title=title, theme="light")
        elif fmt in ("pdf", "png", "svg"):
            gc.plot(res, f"{out}.{fmt}", title=title, theme=args.theme, top=args.top)
    print(f"{len(res.significant)} significant concepts; {len(res.study)} genes used, {len(res.missing)} missing -> {out}.*")


def cmd_compare(args):
    from .compare import draw_compare
    gc = _engine(args)
    bg = _read_ids(args.background)
    sets = _read_sets(args.sets)
    results = {name: gc.run(g, background=bg, **_run_opts(args)) for name, g in sets.items()}
    out = Path(args.out)
    for name, r in results.items():
        stem = f"{out}_{''.join(ch if ch.isalnum() else '_' for ch in name)}"
        r.to_tsv(f"{stem}.tsv")
        gc.plot(r, f"{stem}.pdf", title=f"GeneCloud · {name}", theme=args.theme, top=args.top)
        gc.html(r, f"{stem}.html", title=f"GeneCloud · {name}")
    d = draw_compare(results, f"{out}_compare.pdf")
    draw_compare(results, f"{out}_compare.png")
    d.to_csv(f"{out}_compare.tsv", sep="\t", index=False)
    print(f"{len(results)} sets -> {out}_*")


def cmd_legacy(args):
    from .legacy import legacy_genecloud
    gc = _engine(args)
    t = legacy_genecloud(gc, _read_ids(args.genes), background=_read_ids(args.background),
                         sampling=args.sampling)
    t.to_csv(args.out, sep="\t", index=False)
    print(args.out)


def main(argv=None):
    p = argparse.ArgumentParser(prog="genecloud", description="Semantic enrichment word clouds for gene lists.")
    p.add_argument("--version", action="version", version=f"genecloud {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="download sources and build the annotation table")
    b.add_argument("--species", default="arabidopsis")
    b.add_argument("--sources", default="genecloud_sources", help="folder for the raw source files")
    b.add_argument("--out", default="arabidopsis.jsonl.gz")
    b.add_argument("--offline", action="store_true", help="use files already in --sources")
    b.set_defaults(func=cmd_build)

    def common(sp):
        sp.add_argument("-a", "--annotation", required=True, help="table made by 'genecloud build'")
        sp.add_argument("--species", default="arabidopsis")
        sp.add_argument("-b", "--background", required=True,
                        help="REQUIRED: file with the background genes, i.e. every gene that could have been "
                             "in the list (all expressed genes, all genes on the array...)")
        sp.add_argument("--layers", default="word,phrase,keyword", help="any of word,phrase,keyword,go")
        sp.add_argument("--fdr", type=float, default=0.05)
        sp.add_argument("--min-genes", type=int, default=2)
        sp.add_argument("--max-frac", type=float, default=0.25, help="skip concepts carried by > this fraction of background")
        sp.add_argument("--merge", type=float, default=0.75, help="Jaccard threshold to merge redundant concepts")
        sp.add_argument("--trend-p", type=float, default=0.01)
        sp.add_argument("--theme", default="light", choices=["light", "dark"])
        sp.add_argument("--top", type=int, default=60)

    r = sub.add_parser("run", help="one gene list -> cloud (pdf/png/svg), interactive html and table")
    r.add_argument("genes", help="file with gene ids ('-' for stdin)")
    common(r)
    r.add_argument("-o", "--out", default="genecloud")
    r.add_argument("--title")
    r.add_argument("--formats", default="pdf,png,html")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("compare", help="several gene sets (two-column file gene<TAB>set) -> clouds + dot plot")
    c.add_argument("sets")
    common(c)
    c.add_argument("-o", "--out", default="genecloud")
    c.set_defaults(func=cmd_compare)

    l = sub.add_parser("legacy", help="original 2013 Monte Carlo GeneCloud test on the modern annotation")
    l.add_argument("genes")
    common(l)
    l.add_argument("--sampling", type=int, default=100)
    l.add_argument("-o", "--out", default="genecloud_legacy.tsv")
    l.set_defaults(func=cmd_legacy)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
