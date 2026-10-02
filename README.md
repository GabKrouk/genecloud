# GeneCloud 2

**Semantic enrichment word clouds for gene lists.** Give GeneCloud a list of genes (a cluster, a set of
differentially expressed genes…) and it tells you, in words, what these genes have in common: which words and
phrases of their functional annotation are over-represented compared with the genome, with proper statistics,
drawn as a word cloud in which **the size of each word is its statistical enrichment**.

GeneCloud 2 is a complete rewrite of GeneCloud (G. Krouk, 2013, R package `GeneCloud.gb`) with a modern,
openly licensed genome annotation and a modern statistical engine.

![Auxin-induced genes](examples/images/hormone_IAA_up.png)

## What is new compared with GeneCloud 2013

| | GeneCloud (2013) | GeneCloud 2 |
|---|---|---|
| Annotation | TAIR gene descriptions frozen in Dec 2012, bundled | Rebuilt on demand from **UniProtKB** (protein names, curated FUNCTION text, keywords, families), **NCBI Gene** (symbols, descriptions, nomenclature names) and the **Gene Ontology** (GAF + ontology). All CC BY 4.0 / public domain |
| Unit of counting | every occurrence of a word (a long description weighs more) | **genes** carrying the concept (one gene = one vote) |
| Concepts | single words | words (lemmatised: *transporters* → *transporter*), **two-word phrases** detected genome-wide by normalised PMI (*high affinity*, *plasma membrane*, *xenobiotic detoxification*), **UniProt keywords**, optionally **GO terms** (with ancestors) |
| Statistics | 100 random gene lists drawn from the whole genome; significant if never reached (P < 0.01 at best), no multiple-testing correction | exact **hypergeometric test** against a **mandatory background** (the genes that could have been in the list), **Benjamini–Hochberg FDR**, fold enrichment |
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

# 2. one gene list + its background -> cloud + interactive report + table
genecloud run my_cluster.txt -a arabidopsis.jsonl.gz -b expressed_genes.txt -o results/cluster1

# 3. several lists at once (two-column file: gene <TAB> set) -> one cloud per set + comparison dot plot
genecloud compare clusters.tsv -a arabidopsis.jsonl.gz -b expressed_genes.txt -o results/clusters
```

The background (`-b/--background`, or `background=` in Python) is **required**: it is the list of genes
that *could* have appeared in your list — every gene detected in the RNA-seq, every gene on the array.

Gene files can contain one identifier per line or any text: identifiers are recognised by pattern
(`AT1G01010`, `At1g01010.1` …); lines starting with `#` are ignored.

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

The background is mandatory. Pass the genes that *could* have been in your list: all genes detected in the
RNA-seq, or all genes represented on the microarray. Using the whole genome instead inflates enrichment for
anything associated with expression in your tissue or with presence on the array.

### Good practice

* Small sets (< 15 genes) rarely reach significance after FDR correction; look at the trends, and at the
  `genes` column: a concept carried by 2 genes is a lead, not a conclusion.
* GeneCloud complements, not replaces, GO enrichment: its strength is to surface **vocabulary that no ontology
  term captures** (a protein family, a transported molecule, a phenotype described in UniProt).

## Examples: the lists of the GeneCloud paper

`examples/` re-runs the three analyses of the original GeneCloud paper
(Krouk, Carré, Fizames, Gojon, Ruffel & Lacombe, 2015, *Mol. Plant* 8:971,
[doi:10.1016/j.molp.2015.02.005](https://doi.org/10.1016/j.molp.2015.02.005)) with GeneCloud 2.
All lists come from ATH1 microarrays, so the background is every gene on the ATH1 array
(`examples/lists/background_ATH1.txt`, 21,521 annotated genes).

```bash
genecloud build --species arabidopsis --out arabidopsis.jsonl.gz
python examples/run_examples.py arabidopsis.jsonl.gz                # -> examples/results/, examples/images/
python examples/prepare_lists.py SOURCE_DIR                         # optional: rebuild the lists from the public files
```

| list | source | genes | 2015 paper found | GeneCloud 2 top concepts (FDR) |
|---|---|---|---|---|
| auxin (IAA)-induced | Nemhauser et al. 2006 *Cell*, Table S5 | 430 | *auxin* (P = 3·10⁻¹⁴) | **auxin responsive** (10⁻³⁹), Auxin signaling pathway, promoter element (AuxRE), Aux/IAA, SAUR, acetic acid |
| induced by PHR1 over-expression | Bustos et al. 2010 *PLoS Genet*, GSE20955: OxPHR1 vs *phr1*, ≥ 2-fold, t-test P < 0.05 | 381 | *phosphate*, *purple*, *glutaredoxin*, *phosphatase* | **phosphate starvation** (10⁻²¹), inorganic, purple acid (phosphatase), phosphatase, glycolipid / DGDG / MGDG, monothiol glutaredoxin |
| nitrate-responsive | Canales et al. 2014 *Front. Plant Sci.*, Table S1 (meta-analysis) | 2,264 | *nitrate*, *shaqkyf* (Wang et al. 2004 lists) | hypoxia, photosynthesis, apoplast, **nitrate** (10⁻¹⁸; nitrite, nitrate assimilation, high-affinity nitrate merged), sulfate, ammonium |

![PHR1](examples/images/PHR1_induced.png)

![Nitrate](examples/images/nitrate_responsive.png)

The seven hormone treatments of Nemhauser et al. (2006), induced genes, compared in one dot plot: each hormone
is recognised by its own vocabulary (auxin, brassinosteroid → cell wall / xyloglucan, cytokinin, abscisic acid,
ethylene, jasmonic acid / glucosinolate). Gibberellin (40 genes) gives little, as in the original study.

![Hormones](examples/images/hormones_compare.png)

Notes on the examples

* The 2015 nitrate example used the WT and NR-null lists of Wang et al. (2004, *Plant Physiol.* 136:2512), whose
  supplementary tables are not openly downloadable. It is replaced here by the nitrate-responsive genes of the
  Canales et al. (2014) meta-analysis of 27 ATH1 nitrate-treatment experiments.
* *shaqkyf* came from the 2012 TAIR descriptions of GARP/MYB-related proteins ("myb-like HTH transcriptional
  regulator family protein … SHAQKYF class"). That wording is absent from today's UniProt / NCBI / GO
  annotation, so the term cannot appear any more.
* Auxin is reported for auxin-regulated genes partly because they were *described* from such experiments: the
  circularity discussed in the 2015 paper still applies.
* `examples/results/` holds the tables (`*.tsv`) and interactive reports (`*.html`) of every example.

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

If you use GeneCloud, please cite Krouk G, Carré C, Fizames C, Gojon A, Ruffel S, Lacombe B (2015) GeneCloud
reveals semantic enrichment in lists of gene descriptions. *Mol. Plant* 8:971–973, this repository (see
`CITATION.cff`), and the annotation sources:
UniProt Consortium (Nucleic Acids Res.), Gene Ontology Consortium (Genetics), NCBI Gene (Nucleic Acids Res.).

## License

MIT for the code. Annotation data are downloaded at build time from their providers and keep their own licences
(UniProt and Gene Ontology: CC BY 4.0; NCBI Gene: public domain).
