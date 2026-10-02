"""From annotation text to 'concepts' (words and phrases), one set per gene.

Design choices (vs. GeneCloud 2013):
  * concepts are counted once per gene (presence), so one long description cannot dominate a cloud;
  * words are lower-cased and lightly lemmatised (transporters -> transporter, families -> family);
  * hyphenated and slashed compounds are split ("auxin-responsive" -> auxin, responsive);
  * identifiers, numbers, PubMed / ECO evidence tags and annotation boilerplate are removed;
  * recurrent two-word phrases ("nitrate transport", "cell wall") are detected genome-wide by normalised
    pointwise mutual information and become concepts of their own;
  * UniProt keywords and GO terms can be added as whole concepts (layers).
"""
from __future__ import annotations

import math
import re
from collections import Counter
from typing import Iterable

ENGLISH_STOPWORDS = frozenset("""
a about above after again against all also although am among an and any are as at be because been before being
below between both but by can cannot could did do does doing down during each either else etc even ever every few
for from further had has have having he her here hers herself him himself his how however i if in into is it its
itself just least less many may me might more most much must my myself neither no nor not now of off often on once
only or other others otherwise our ours ourselves out over own per rather same several shall she should since so
some such than that the their theirs them themselves then there therefore these they this those though through
thus to too under until up upon us very via was we were what when where whether which while who whom whose why
will with within without would yet you your yours yourself yourselves
""".split())

# words that occur everywhere in gene annotations and carry no biological signal
ANNOTATION_STOPWORDS = frozenset("""
protein proteins gene genes putative probable possible possibly likely uncharacterized unknown hypothetical
family families member members subfamily superfamily like related similar similarity homolog homologue homolog
homology ortholog orthologue paralog domain domains containing contains contain contained motif region type
class group subunit activity activities process processes involved involve involves involving function functions
functional functioning role roles play plays playing required requires require acts act acting mediate mediates
mediated mediating encode encodes encoded encoding expressed expression express expresses level levels
specific specifically associated associates association interact interacts interacting interaction
interactions bind binds binding regulated show shows shown display displays may also via well due including
include includes thereby part essential important major minor main several various different number numerous
form forms formation present presence product products term terms cellular cell-specific evidence predicted
prediction sequence sequences reference proteome isoform isoforms fragment complete partial chain chains
probably required unclear least other cause causes caused result results resulting lead leads leading
positively negatively positive negative component components factor factors
previously recently respectively directly indirectly together mainly mostly partly partially highly strongly
weakly especially particularly additionally subsequently generally usually normally whereas across another
along around based known described identified suggested thought believed found named called termed designated
first second third one two three four five many multiple single several either upon
mutant mutants mutation mutations allele alleles phenotype phenotypes wild wild-type plant plants line lines
transgenic knockout overexpression overexpressing overexpressor observed compared comparison exhibit exhibits
exhibited identical encoding source study studies data analysed analyzed confirmed demonstrated indicate
indicates indicated suggest suggests reported revealed journal
""".split())
# note: 'regulation', 'transport', 'response', 'binding' are kept on purpose: in phrases they are informative.
ANNOTATION_STOPWORDS = ANNOTATION_STOPWORDS - {"binding"}

_EVIDENCE = re.compile(r"\{ECO:[^}]*\}|\((?:PubMed|By similarity)[^)]*\)|PubMed:\d+", re.I)
_EC = re.compile(r"\bEC\s?\d+(\.[\d-]+){0,3}\b")
_TOKEN = re.compile(r"[a-z][a-z0-9'\u2011]*[a-z0-9]|[a-z]")
_PREFIX = re.compile(r"\b(trans|cis|co|non|anti|multi|pre|post|pro|sub|semi|mono|poly|inter|intra|extra|hetero|homo|self|cross|short)-(?=[a-z])")
_HAS_DIGIT = re.compile(r"\d")
_LEMMA_KEEP = frozenset("""species series analysis biosynthesis synthesis hydrolysis homeostasis apoptosis
mitosis meiosis stress process class glass gas lens chitosanase mucus virus focus status bus axis basis
thesis genesis osmosis lysis pathogenesis embryogenesis morphogenesis organogenesis senescence nucleus
cactus pius plus thus has was is its this tris abscisic trans""".split())


def lemma(w: str) -> str:
    """Very light English lemmatiser for annotation vocabulary (plural -> singular)."""
    if len(w) <= 4 or w in _LEMMA_KEEP or w.endswith(("sis", "ss", "us", "ous", "ics")):
        return w
    if w.endswith("ies") and len(w) > 5:
        return w[:-3] + "y"
    if w.endswith(("ches", "shes", "xes", "sses")):
        return w[:-2]
    if w.endswith("s") and not w.endswith(("is", "as")):
        return w[:-1]
    return w


class Tokenizer:
    def __init__(self, id_pattern: re.Pattern | None = None, extra_stopwords: Iterable[str] = (),
                 keep_symbols: bool = False, min_len: int = 3):
        self.id_pattern = id_pattern
        self.stop = ENGLISH_STOPWORDS | ANNOTATION_STOPWORDS | {w.lower() for w in extra_stopwords}
        self.keep_symbols = keep_symbols
        self.min_len = min_len

    def sentences(self, text: str) -> list[list[str]]:
        """Token runs, broken at punctuation and stop words (used for phrase detection)."""
        if not text:
            return []
        text = _EVIDENCE.sub(" ", text)
        text = _EC.sub(" ", text)
        if self.id_pattern is not None:
            text = self.id_pattern.sub(" ", text)
        text = _PREFIX.sub("\\1\u2011", text.lower())          # keep 'trans-golgi', 'cis-acting' as one token
        text = text.replace("-", " ").replace("/", " ")
        runs, cur = [], []
        for chunk in re.split(r"[.,;:()\[\]{}\"]+", text):
            for raw in chunk.split():
                m = _TOKEN.fullmatch(raw.strip("'"))
                w = m.group(0) if m else None
                ok = (w is not None and len(w) >= self.min_len and w not in self.stop
                      and (self.keep_symbols or not _HAS_DIGIT.search(w)))
                if ok:
                    cur.append(lemma(w))
                elif cur:
                    runs.append(cur)
                    cur = []
            if cur:
                runs.append(cur)
                cur = []
        return runs

    def words(self, text: str) -> set[str]:
        return {w for run in self.sentences(text) for w in run}


def detect_phrases(gene_runs: dict[str, list[list[str]]], min_genes: int = 5, min_npmi: float = 0.35,
                   max_phrases: int = 20000) -> set[tuple[str, str]]:
    """Two-word collocations, counted at gene level, scored by normalised PMI."""
    uni, bi = Counter(), Counter()
    for runs in gene_runs.values():
        u, b = set(), set()
        for r in runs:
            u.update(r)
            b.update(zip(r, r[1:]))
        uni.update(u)
        bi.update(b)
    n = max(len(gene_runs), 1)
    scored = []
    for (a, b), c in bi.items():
        if c < min_genes or a == b:
            continue
        p_ab, p_a, p_b = c / n, uni[a] / n, uni[b] / n
        npmi = math.log(p_ab / (p_a * p_b)) / -math.log(p_ab)
        if npmi >= min_npmi:
            scored.append((npmi, (a, b)))
    scored.sort(reverse=True)
    return {p for _, p in scored[:max_phrases]}


def phrases_in(runs: list[list[str]], phrases: set[tuple[str, str]]) -> set[str]:
    return {f"{a} {b}" for r in runs for a, b in zip(r, r[1:]) if (a, b) in phrases}
