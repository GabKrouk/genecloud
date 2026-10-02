# GeneCloud 2

**Semantic enrichment word clouds for gene lists.** Give GeneCloud a list of genes (a cluster, a set of
differentially expressed genes…) and it tells you, in words, what these genes have in common: which words and
phrases of their functional annotation are over-represented compared with the genome, with proper statistics,
drawn as a word cloud in which **the size of each word is its statistical enrichment**.

GeneCloud 2 is a complete rewrite of GeneCloud (G. Krouk, 2013, R package `GeneCloud.gb`) with a modern,
openly licensed genome annotation and a modern statistical engine.



## What is new compared with GeneCloud 2013

| | GeneCloud (2013) | GeneCloud 2 |
|---|---|---|
| Annotation | TAIR gene descriptions frozen in Dec 2012, bundled | Rebuilt on demand from **UniProtKB** (protein names, curated FUNCTION text, keywords, families), **NCBI Gene** (symbols, descriptions, nomenclature names) and the **Gene Ontology** (GAF + ontology). All CC BY 4.0 / public domain |
| Unit of counting | every occurrence of a word (a long description weighs more) | **genes** carrying the concept (one gene = one vote) |
| Concepts | single words | words (lemmatised: *transporters* → *transporter*), **two-word phrases** detected genome-wide by normalised PMI (*high affinity*, *plasma membrane*, *xenobiotic detoxification*), **UniProt keywords**, optionally **GO terms** (with ancestors) |
| Statistics | 100 random gene lists; significant if never reached (P < 0.01 at best), no multiple-testing correction | exact **hypergeometric test** against a user-defined **background** (e.g. all expressed genes), **Benjamini–Hochberg FDR**, fold enrichment |
| Noise | gene names, numbers and identifiers appear as "significant" | identifiers, numbers, evidence tags and annotation boilerplate removed; a concept must be carried by ≥ 2 genes; over-general concepts (> 25 % of background) are not tested |
| Redundancy | *ammonium*, *ammonium transporter*, *amt11* … all shown | redundant concepts carried by the same genes are **merged** (the most significant represents the group); lone words are shown with the phrase their genes share (*distance* → *long distance*) |
| Word size | word frequency | **−log10 FDR**; colour = fold enrichment; non-significant trends optionally in grey |
| Layout | random | deterministic spiral layout (same input → same picture) |
| Outputs | PDF | PDF / PNG / SVG, a **self-contained interactive HTML** report (hover a word for its statistics, click it to see its genes), TSV tables, and a **multi-set comparison** dot plot |
| Reproducibility | — | the original 2013 test is kept (`genecloud legacy`) to compare both methods on the same annotation |

## Install

```bash
pip install git+https://github.com/GabKrouk/genecloud
# or, from a clone:
pip install -e ".[test]"
```

Python ≥ 3.9; dependencies: numpy, pandas, scipy, matplotlib.

## Quick start

```bash
# 1. build the annotation once (~30 MB download, ~10 s to parse)
genecloud build --species arabidopsis --out arabidopsis.jsonl.gz

# 2. one gene list -> cloud + interactive report + table
genecloud run my_cluster.txt -a arabidopsis.jsonl.gz --background expressed_genes.txt -o results/cluster1

# 3. several lists at once (two-column file: gene <TAB> set) -> one cloud per set + comparison dot plot
genecloud compare clusters.tsv -a arabidopsis.jsonl.gz --background expressed_genes.txt -o results/clusters
```

Gene files can contain one identifier per line or any text: identifiers are recognised by pattern
(`AT1G01010`, `At1g01010.1` …).

From Python:

```python
from genecloud import GeneCloud

gc = GeneCloud("arabidopsis.jsonl.gz", layers=("word", "phrase", "keyword"))
res = gc.run(genes, background=expressed, fdr=0.05, min_genes=2)
res.significant            # pandas DataFrame: concept, k/n, K/N, fold, P, FDR, genes, merged synonyms
gc.plot(res, "cloud.pdf", title="Cluster 2")
gc.html(res, "cloud.html", title="Cluster 2")
```

## How it works

1. **Annotation.** For every gene, GeneCloud gathers its NCBI description and nomenclature names, UniProt
   protein names, curated FUNCTION paragraph, families and keywords, and the names of its direct GO annotations.
   Evidence tags (`PubMed:…`, `{ECO:…}`), EC numbers and gene identifiers are stripped.
2. **Concepts.** Text is lower-cased, split at punctuation and stop words (English + annotation boilerplate such
   as *protein*, *putative*, *involved*, *family*), and lightly lemmatised. Compound prefixes are kept
   (*trans-Golgi*). Two-word phrases that co-occur far more than expected across the genome (normalised PMI ≥ 0.35,
   ≥ 5 genes) become concepts too. Each gene is reduced to a **set** of concepts.
3. **Test.** For each concept carried by *k* of the *n* study genes and *K* of the *N* background genes,
   P = hypergeometric upper tail; FDR by Benjamini–Hochberg over all concepts carried by ≥ `min_genes` study genes.
4. **Redundancy.** Significant concepts are visited from the most to the least significant; a concept joins an
   earlier one when their study genes overlap with Jaccard ≥ 0.75, or when the words of one are contained in the
   other and they share most genes. The cloud shows one word per group; the table lists the merged synonyms.
5. **Cloud.** Word size ∝ −log10 FDR, colour = fold enrichment (colour-blind-friendly ramp), italic = UniProt
   keyword, grey = trend (P ≤ 0.01 but FDR above the cut-off). Words are placed on a deterministic spiral and the
   font scale shrinks until every significant concept fits.

### Choosing the background

Always pass the genes that *could* have been in your list (e.g. all genes detected in the RNA-seq). Using the
whole genome as background inflates enrichment for anything associated with expression in your tissue.

### Good practice

* Small sets (< 15 genes) rarely reach significance after FDR correction; look at the trends, and at the
  `genes` column: a concept carried by 2 genes is a lead, not a conclusion.
* GeneCloud complements, not replaces, GO enrichment: its strength is to surface **vocabulary that no ontology
  term captures** (a protein family, a transported molecule, a phenotype described in UniProt).

## Other species

Species are declared in `src/genecloud/species.py` (GAF URL, NCBI gene_info URL, UniProt reference proteome,
identifier pattern). Adding one is a few lines; pull requests welcome.

## Output columns (`*.tsv`)

| column | meaning |
|---|---|
| `layer` | word, phrase, keyword or GO |
| `label` / `display` | the concept / the text shown in the cloud |
| `k`, `n`, `K`, `N` | study genes with the concept, study size, background genes with the concept, background size |
| `fold` | (k/n)/(K/N) |
| `p`, `fdr` | hypergeometric P, Benjamini–Hochberg FDR |
| `status` | significant / trend / (empty) |
| `representative`, `group`, `merged` | redundancy grouping |
| `genes`, `symbols` | study genes carrying the concept |

## Citing

If you use GeneCloud, please cite this repository (see `CITATION.cff`) and the annotation sources:
UniProt Consortium (Nucleic Acids Res.), Gene Ontology Consortium (Genetics), NCBI Gene (Nucleic Acids Res.).

## License

MIT for the code. Annotation data are downloaded at build time from their providers and keep their own licences
(UniProt and Gene Ontology: CC BY 4.0; NCBI Gene: public domain).
