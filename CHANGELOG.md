# Changelog

## 2.0.0 — 2026-10-02
Complete rewrite in Python of GeneCloud (R package `GeneCloud.gb`, 2013).

- Annotation rebuilt on demand from UniProtKB, NCBI Gene and the Gene Ontology (open licences).
- Gene-level counting, lemmatisation, genome-wide phrase detection (normalised PMI), UniProt keywords and GO layers.
- Hypergeometric test against a custom background, Benjamini–Hochberg FDR, fold enrichment.
- Redundant concepts merged; lone words displayed with their shared phrase.
- Word size = −log10 FDR, colour = fold enrichment, deterministic layout, light and dark themes.
- Interactive self-contained HTML report, TSV tables, multi-set comparison dot plot.
- `genecloud legacy` reproduces the 2013 Monte Carlo test on the modern annotation.
