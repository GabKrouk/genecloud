"""Rebuild the example gene lists of Krouk et al. (2015, Mol. Plant 8:971) from the public source files.

Source files (download them into one folder and pass it as argument):

  nem_mmc2.xls ... nem_mmc8.xls   Nemhauser, Hong & Chory (2006) Cell 126:467, Supplementary Tables 1-7
                                  https://doi.org/10.1016/j.cell.2006.05.050 (Supplemental Data, mmc2-mmc8)
  GSE20955_series_matrix.txt.gz   Bustos et al. (2010) PLoS Genet 6:e1001102, OxPHR1 vs phr1 (-Pi, DEX, CHX)
                                  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20955/matrix/
  GPL198.annot.gz                 Affymetrix ATH1 annotation (probe -> AGI)
                                  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL0nnn/GPL198/annot/
  canales_S1.xls                  Canales et al. (2014) Front Plant Sci 5:22, Table S1 (Data Sheet 1)
                                  https://doi.org/10.3389/fpls.2014.00022

.xls files are read with pandas (needs `xlrd`) or, failing that, converted with LibreOffice.

Usage:  python examples/prepare_lists.py SOURCE_DIR  [OUT_DIR=examples/lists]
"""
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

AGI = re.compile(r"^AT[1-5CM]G\d{5}$")
HORMONES = {"mmc2": "ABA", "mmc3": "ACC", "mmc4": "BL", "mmc5": "GA", "mmc6": "IAA", "mmc7": "MeJA", "mmc8": "CK"}


def read_xls(path, **kw):
    try:
        return pd.read_excel(path, header=None, **kw)
    except ImportError:
        tmp = Path(tempfile.mkdtemp())
        subprocess.run(["soffice", "--headless", "--convert-to", "xlsx", "--outdir", str(tmp), str(path)],
                       check=True, capture_output=True)
        return pd.read_excel(tmp / (Path(path).stem + ".xlsx"), header=None, **kw)


def agis(values):
    out = []
    for v in values:
        v = str(v).strip().upper()
        if AGI.match(v) and v not in out:
            out.append(v)
    return out


def ath1_map(annot):
    with gzip.open(annot, "rt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("#", "!", "^"))]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str)
    orf = df["Platform_ORF"].fillna("").str.upper().str.strip()
    # keep probes that match exactly one AGI
    ok = orf.str.match(AGI.pattern)
    return dict(zip(df.loc[ok, "ID"], orf[ok]))


def nemhauser(src, out):
    for mmc, hormone in HORMONES.items():
        t = read_xls(src / f"nem_{mmc}.xls")
        hdr = t.index[t.iloc[:, 0].astype(str).str.lower().eq("up")][0]
        body = t.iloc[hdr + 1:]
        cols = {str(v).lower(): j for j, v in enumerate(t.iloc[hdr]) if str(v).lower() in ("up", "down", "complex")}
        for direction in ("up", "down"):
            g = agis(body.iloc[:, cols[direction]])
            write(out / f"hormone_{hormone}_{direction}.txt", g,
                  f"Nemhauser et al. 2006, {hormone}-regulated genes ({direction}), linear model P<0.01")


def phr1(src, out, probe2agi):
    with gzip.open(src / "GSE20955_series_matrix.txt.gz", "rt") as fh:
        lines = fh.read().split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    from io import StringIO
    m = pd.read_csv(StringIO(lines), sep="\t", index_col=0)
    ox, mut = m.iloc[:, :3], m.iloc[:, 3:6]           # GSM523974-76 OxPHR1, GSM523977-79 phr1 (log2 RMA)
    lfc = ox.mean(axis=1) - mut.mean(axis=1)
    p = stats.ttest_ind(ox, mut, axis=1).pvalue
    sel = m.index[(lfc >= 1) & (p < 0.05)]
    g = agis(probe2agi.get(pr, "") for pr in sel)
    write(out / "PHR1_induced.txt", g,
          "Bustos et al. 2010 (GSE20955): >=2-fold induced in OxPHR1 vs phr1 (-Pi, DEX+CHX), t-test P<0.05")


def nitrate(src, out):
    t = read_xls(src / "canales_S1.xls")
    g = agis(t.iloc[1:, 1])
    write(out / "nitrate_responsive.txt", g,
          "Canales et al. 2014, Table S1: nitrate-responsive genes across 27 ATH1 experiments (meta-analysis)")


def write(path, genes, comment):
    path.write_text(f"# {comment}\n# {len(genes)} genes\n" + "\n".join(genes) + "\n")
    print(f"{path.name:32s} {len(genes):5d}")


def main(argv=None):
    argv = argv or sys.argv[1:]
    src = Path(argv[0])
    out = Path(argv[1]) if len(argv) > 1 else Path(__file__).parent / "lists"
    out.mkdir(parents=True, exist_ok=True)
    probe2agi = ath1_map(src / "GPL198.annot.gz")
    bg = sorted(set(probe2agi.values()))
    write(out / "background_ATH1.txt", bg, "All AGI represented by a single-gene probe set on the Affymetrix ATH1 array (GPL198)")
    nemhauser(src, out)
    phr1(src, out, probe2agi)
    nitrate(src, out)


if __name__ == "__main__":
    main()
