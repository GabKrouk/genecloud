# Changelog

## 2.1.0 — 2026-10-02 · "GeneCloud annotation 2026"
- **Annotation 2026**: the latest TAIR public data release (TAIR_Data_20250930: Araport11 short descriptions,
  curator summaries, computational descriptions, gene aliases, germplasm and allele phenotypes, Araport11 gene
  types) is downloaded and merged with UniProtKB, NCBI Gene and GO. Text presets `all` (default), `tair`, `open`.
- **Web application** (`web/`, GitHub Pages): the 2015 GeneCloud page re-done as a static site that runs the
  analysis in the browser; `genecloud export-web` writes its data. The JavaScript engine reproduces the Python
  results exactly (tested).
- Statistics: hypothesis family fixed by the background (min_genes ≤ K ≤ 25 % N, absent concepts with P = 1);
  concepts below min_genes enter with P = 1; Benjamini–Yekutieli option (`--adjust BY`); P and FDR computed in log
  scale (`log10_p`, `log10_fdr`), exact below 1e-300; deterministic ordering.
- Redundancy: groups keyed by concept (a phrase and a GO term with the same name are no longer mixed); nested
  wording merges only when half of the larger gene set is shared; a lone word takes a phrase label only when
  ≥ 90 % of its genes carry the phrase.
- Cloud: size key and colour key, `--size-by fold`; fixed words overlapping the subtitle and the colour-bar label
  falling outside the figure; `scripts/check_layout.py`.
- Result objects report `not_in_background` and `not_annotated` separately. Lemmatiser fix (*catalyzes*), more
  annotation boilerplate stop words (*mutant*, *journal*, …).

## 2.0.1 — 2026-10-02
- The background gene list is now **mandatory** (`-b/--background` in the CLI, `background=` in `GeneCloud.run`
  and `legacy_genecloud`).
- `examples/`: the gene lists of the 2015 GeneCloud paper (Nemhauser 2006 hormones, PHR1 / GSE20955, nitrate)
  rebuilt from the public data, an ATH1 background, the scripts and the resulting clouds, tables and reports.
- Lines starting with `#` are ignored in gene files; `short-` compounds kept as one token; concept cache
  invalidated when the text-processing code changes.

## 2.0.0 — 2026-10-02
Complete rewrite in Python of GeneCloud (R package `GeneCloud.gb`, 2013).

- Annotation rebuilt on demand from UniProtKB, NCBI Gene and the Gene Ontology (open licences).
- Gene-level counting, lemmatisation, genome-wide phrase detection (normalised PMI), UniProt keywords and GO layers.
- Hypergeometric test against a custom background, Benjamini–Hochberg FDR, fold enrichment.
- Redundant concepts merged; lone words displayed with their shared phrase.
- Word size = −log10 FDR, colour = fold enrichment, deterministic layout, light and dark themes.
- Interactive self-contained HTML report, TSV tables, multi-set comparison dot plot.
- `genecloud legacy` reproduces the 2013 Monte Carlo test on the modern annotation.
