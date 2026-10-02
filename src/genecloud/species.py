"""Species presets: where to download the annotation sources and how to recognise gene identifiers.

Every source is openly licensed:
  * Gene Ontology annotations (GAF) and ontology  - CC BY 4.0, http://geneontology.org
  * UniProtKB reference proteome                  - CC BY 4.0, https://www.uniprot.org
  * NCBI Gene (gene_info)                         - public domain, https://www.ncbi.nlm.nih.gov/gene
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Species:
    name: str
    taxon: int
    gaf_url: str
    gene_info_url: str
    uniprot_proteome: str
    id_regex: str                      # matches a locus identifier anywhere in a string
    id_case: str = "upper"             # how identifiers are normalised
    extra_stopwords: tuple = field(default_factory=tuple)

    @property
    def id_pattern(self) -> re.Pattern:
        return re.compile(self.id_regex, re.IGNORECASE)

    def normalise(self, gene_id: str) -> str:
        g = gene_id.strip().split(".")[0]          # drop isoform suffix (AT1G01010.1)
        return g.upper() if self.id_case == "upper" else g

    @property
    def uniprot_url(self) -> str:
        fields = ",".join(["accession", "reviewed", "gene_oln", "gene_primary", "gene_synonym",
                           "protein_name", "cc_function", "keyword", "protein_families"])
        return ("https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv"
                f"&query=%28proteome%3A{self.uniprot_proteome}%29&fields={fields}")


GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"

SPECIES = {
    "arabidopsis": Species(
        name="Arabidopsis thaliana",
        taxon=3702,
        gaf_url="https://current.geneontology.org/annotations/tair.gaf.gz",
        gene_info_url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/Arabidopsis_thaliana.gene_info.gz",
        uniprot_proteome="UP000006548",
        id_regex=r"\bAT[1-5CM]G\d{5}\b",
        extra_stopwords=("arabidopsis", "thaliana", "tair"),
    ),
}
