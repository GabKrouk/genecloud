// GeneCloud 2 — statistical engine, JavaScript port of genecloud.core.GeneCloud.run (Python).
// Runs in the browser and in Node (tests/test_web_engine.mjs checks it against the Python results).
//
//   universe  = genes annotated in the selected layers ∩ background            (N genes)
//   family    = concepts with min_genes <= K <= max_concept_frac * N            (m hypotheses, background only)
//   P         = hypergeometric upper tail P(X >= k), X ~ H(N, K, n)
//   FDR       = Benjamini–Hochberg (or –Yekutieli) over the m hypotheses, in log10 (no underflow); concepts
//               absent from the list, or carried by fewer than min_genes genes, enter with P = 1
//   status    = significant (k >= min_genes, FDR <= fdr, fold >= min_fold) | trend (P <= trend_p)
//   redundancy: greedy merging of concepts carried by the same genes (Jaccard) or nested wording.

const LAYER_NAMES = { w: "word", p: "phrase", k: "keyword", b: "GO BP", f: "GO MF", c: "GO CC", g: "GO" };
const LAYER_GROUP = { w: "word", p: "phrase", k: "keyword", b: "go", f: "go", c: "go", g: "go" };
const LAYER_BIT = { word: 1, phrase: 2, keyword: 4, go: 8 };
const NBH = "‑";

// ------------------------------------------------------------------------------------------- loading
export function decodeIndex(bytes, nConcepts) {
  // varint(count), then delta-encoded gene indices, for each concept
  const start = new Uint32Array(nConcepts + 1);
  let total = 0, pos = 0;
  const readVar = () => { let x = 0, s = 0, b; do { b = bytes[pos++]; x += (b & 127) * 2 ** s; s += 7; } while (b & 128); return x; };
  // first pass: counts
  const counts = new Uint32Array(nConcepts);
  const save = [];
  for (let c = 0; c < nConcepts; c++) {
    const n = readVar(); counts[c] = n; save.push(pos);
    for (let i = 0; i < n; i++) readVar();
    total += n;
  }
  const mem = new Uint32Array(total);
  let o = 0;
  for (let c = 0; c < nConcepts; c++) {
    start[c] = o; pos = save[c]; let prev = 0;
    for (let i = 0; i < counts[c]; i++) { prev += readVar(); mem[o++] = prev; }
  }
  start[nConcepts] = o;
  return { start, mem };
}

export class Annotation {
  /** genes: {ids, symbols, desc, flags, bits:{preset:[...]}}, blocks: {name: {labels, layers, ids?, start, mem}} */
  constructor(meta, genes, blocks) {
    this.meta = meta; this.genes = genes; this.blocks = blocks;
    this.nGenes = genes.ids.length;
    this.index = new Map(genes.ids.map((g, i) => [g, i]));
    this.idRe = new RegExp(meta.id_regex.replace(/^\\b|\\b$/g, ""), "gi");
  }
  /** identifiers in a text, normalised (AT1G01010.1 -> AT1G01010), unique, in order */
  parseIds(text) {
    const out = [], seen = new Set();
    for (const m of String(text).matchAll(new RegExp(`\\b${this.idRe.source}(?:\\.\\d+)?\\b`, "gi"))) {
      const g = m[0].split(".")[0].toUpperCase();
      if (!seen.has(g)) { seen.add(g); out.push(g); }
    }
    return out;
  }
  /** concept table for a text preset + layer selection (cached) */
  concepts(preset, layers) {
    const key = preset + "|" + [...layers].sort().join(",");
    this._cache ??= new Map();
    if (this._cache.has(key)) return this._cache.get(key);
    const parts = [];
    const add = (name, filter) => {
      const b = this.blocks[name];
      for (let i = 0; i < b.labels.length; i++) {
        const code = b.layers[i];
        if (!filter(code)) continue;
        const prefix = code === "w" ? "w:" : code === "p" ? "p:" : code === "k" ? "k:" : "g:";
        parts.push({ block: b, i, code, label: b.labels[i], key: prefix + (b.ids ? b.ids[i] : b.labels[i]),
          start: b.start[i], end: b.start[i + 1] });
      }
    };
    if (layers.has("word") || layers.has("phrase"))
      add("wp_" + preset, c => (c === "w" && layers.has("word")) || (c === "p" && layers.has("phrase")));
    if (layers.has("keyword")) add("kw", () => true);
    if (layers.has("go")) add("go", () => true);
    // forward index gene -> phrase concepts (for the display of lone words)
    const phraseOf = new Map();
    parts.forEach((c, j) => {
      if (c.code !== "p") return;
      for (let x = c.start; x < c.end; x++) {
        const g = c.block.mem[x];
        let a = phraseOf.get(g); if (!a) phraseOf.set(g, a = []); a.push(j);
      }
    });
    const res = { list: parts, phraseOf };
    this._cache.set(key, res);
    return res;
  }
}

// ------------------------------------------------------------------------------------ hypergeometric
let LF = new Float64Array(1);
function logFactorials(n) {
  if (LF.length > n) return LF;
  const t = new Float64Array(n + 1);
  // exact summation of log(i) with Kahan compensation
  let s = 0, comp = 0;
  for (let i = 2; i <= n; i++) { const y = Math.log(i) - comp; const z = s + y; comp = (z - s) - y; s = z; t[i] = s; }
  LF = t; return t;
}
const lchoose = (n, k) => LF[n] - LF[k] - LF[n - k];

/** P(X >= k), X ~ hypergeometric(N population, K successes, n draws). */
export function hypergeomSF(k, N, K, n) {
  logFactorials(N);
  const lo = Math.max(0, n - (N - K)), hi = Math.min(K, n);
  if (k <= lo) return 1;
  if (k > hi) return 0;
  const logpmf = i => lchoose(K, i) + lchoose(N - K, n - i) - lchoose(N, n);
  const mean = n * K / N;
  if (k > mean) {                    // upper tail, terms decrease: sum forward from k
    let t = 1, s = 1;
    for (let i = k; i < hi; i++) {
      t *= ((K - i) * (n - i)) / ((i + 1) * (N - K - n + i + 1));
      s += t;
      if (t < s * 1e-17) break;
    }
    return Math.min(1, Math.exp(logpmf(k) + Math.log(s)));
  }
  // otherwise 1 - P(X <= k-1), summing the lower tail backward from k-1
  let t = 1, s = 1;
  for (let i = k - 1; i > lo; i--) {
    t *= (i * (N - K - n + i)) / ((K - i + 1) * (n - i + 1));
    s += t;
    if (t < s * 1e-17) break;
  }
  return Math.max(0, Math.min(1, 1 - Math.exp(logpmf(k - 1) + Math.log(s))));
}

/** log10 P(X >= k); exact far below the double-precision range. */
export function log10HypergeomSF(k, N, K, n) {
  const p = hypergeomSF(k, N, K, n);
  if (p > 1e-280) return Math.log10(p);
  if (p === 0 && k > Math.min(K, n)) return -Infinity;
  logFactorials(N);
  // tiny P: upper tail, terms decrease from k (k is far above the mean)
  let t = 1, s = 1;
  for (let i = k; i < Math.min(K, n); i++) {
    t *= ((K - i) * (n - i)) / ((i + 1) * (N - K - n + i + 1));
    s += t;
    if (t < s * 1e-17) break;
  }
  return (lchoose(K, k) + lchoose(N - K, n - k) - lchoose(N, n) + Math.log(s)) / Math.LN10;
}

/** Benjamini–Hochberg (or –Yekutieli) in log10; m = total number of hypotheses (missing ones have P = 1). */
export function log10Adjust(lp, m, method = "BH") {
  const n = lp.length;
  m = Math.max(m ?? n, n);
  let c = 0;
  if (method === "BY") { let h = 0; for (let i = 1; i <= m; i++) h += 1 / i; c = Math.log10(h); }
  const o = [...lp.keys()].sort((a, b) => lp[a] - lp[b] || a - b);
  const q = new Float64Array(n);
  let run = Infinity;
  for (let r = n - 1; r >= 0; r--) {
    run = Math.min(run, lp[o[r]] + Math.log10(m) - Math.log10(r + 1) + c);
    q[o[r]] = Math.min(run, 0);
  }
  return q;
}

/** Benjamini–Hochberg; m = total number of hypotheses (missing ones have P = 1). */
export function bh(p, m) {
  const n = p.length;
  m = Math.max(m ?? n, n);
  const o = [...p.keys()].sort((a, b) => p[a] - p[b] || a - b);
  const q = new Float64Array(n);
  let run = Infinity;
  for (let r = n - 1; r >= 0; r--) {
    const v = p[o[r]] * m / (r + 1);
    run = Math.min(run, v);
    q[o[r]] = Math.min(run, 1);
  }
  return q;
}

function stem(w) {
  for (const suf of ["ers", "er", "ions", "ion", "ing", "ed"])
    if (w.endsWith(suf) && w.length - suf.length >= 4) return w.slice(0, -suf.length);
  return w;
}
const pyCmp = (a, b) => (a < b ? -1 : a > b ? 1 : 0);   // code-point order like Python str comparison (BMP)

// --------------------------------------------------------------------------------------------- run
export const DEFAULTS = { preset: "all", layers: ["word", "phrase", "keyword"], minGenes: 2, fdr: 0.05,
  maxFrac: 0.25, mergeJaccard: 0.75, minFold: 1.0, trendP: 0.01, adjust: "BH" };

/**
 * genes, background: arrays of identifiers (background REQUIRED).
 * Returns {rows, study, missing, N, n, m, params}; rows sorted like the Python table.
 */
export function run(ann, genes, background, opts = {}) {
  const o = { ...DEFAULTS, ...opts };
  if (!background || !background.length) throw new Error("A background gene list is required.");
  const layers = new Set(o.layers);
  const mask = [...layers].reduce((a, l) => a | LAYER_BIT[l], 0);
  const bits = ann.genes.bits[o.preset];
  const inU = new Uint8Array(ann.nGenes);
  const bgIds = new Set(background.map(g => g.split(".")[0].toUpperCase()));
  let N = 0;
  for (const g of bgIds) {
    const i = ann.index.get(g);
    if (i !== undefined && (bits[i] & mask) && !inU[i]) { inU[i] = 1; N++; }
  }
  if (!N) throw new Error("None of the background genes is annotated.");
  const req = [...new Set(genes.map(g => g.split(".")[0].toUpperCase()))];
  const study = [], missing = [], notInBackground = [], notAnnotated = [];
  const inS = new Uint8Array(ann.nGenes);
  for (const g of req) {
    const i = ann.index.get(g);
    if (i !== undefined && inU[i]) { study.push(g); inS[i] = 1; continue; }
    missing.push(g);
    (bgIds.has(g) ? notAnnotated : notInBackground).push(g);
  }
  const n = study.length;
  if (!n) throw new Error("None of the genes is annotated / in the background.");
  const { list, phraseOf } = ann.concepts(o.preset, layers);
  const rows = [];
  let m = 0;
  for (let j = 0; j < list.length; j++) {
    const c = list[j], mem = c.block.mem;
    let K = 0, k = 0;
    for (let x = c.start; x < c.end; x++) { const g = mem[x]; if (inU[g]) { K++; if (inS[g]) k++; } }
    if (K < o.minGenes || K > o.maxFrac * N) continue;
    m++;
    if (!k) continue;
    const gs = [];
    for (let x = c.start; x < c.end; x++) if (inS[mem[x]]) gs.push(mem[x]);
    rows.push({ j, concept: c.key, layer: LAYER_NAMES[c.code], label: c.label, k, n, K, N,
      fold: (k / n) / (K / N), log10_p: log10HypergeomSF(k, N, K, n), genes: gs });
  }
  rows.forEach(r => { r.tested = r.k >= o.minGenes; });
  // concepts below min_genes can never be called: they enter the correction with P = 1
  const q = log10Adjust(rows.map(r => r.tested ? r.log10_p : 0), m, o.adjust);
  rows.forEach((r, i) => { r.log10_fdr = q[i]; r.p = 10 ** r.log10_p; r.fdr = 10 ** q[i]; });
  rows.sort((a, b) => a.log10_fdr - b.log10_fdr || a.log10_p - b.log10_p || b.k - a.k || pyCmp(a.concept, b.concept));
  for (const r of rows) {
    r.status = r.tested && r.fdr <= o.fdr && r.fold >= o.minFold ? "significant"
      : r.tested && r.p <= o.trendP && r.fold >= o.minFold ? "trend" : "";
    r.representative = false; r.group = ""; r.display = r.label;
  }
  // redundancy
  const cand = rows.filter(r => r.status);
  const reps = [];
  for (const r of cand) {
    const gs = new Set(r.genes), toks = new Set(r.label.toLowerCase().split(/\s+/).filter(Boolean).map(stem));
    let joined = false;
    for (const [s, rs, rt] of reps) {
      let inter = 0; for (const g of gs) if (rs.has(g)) inter++;
      const jac = inter / (gs.size + rs.size - inter);
      const contained = [...toks].every(t => rt.has(t)) || [...rt].every(t => toks.has(t));
      if ((jac >= o.mergeJaccard || (contained && inter / Math.max(gs.size, rs.size) >= 0.5))
          && (r.status === s.status || s.status === "significant")) {
        r.group = s.concept; joined = true; break;
      }
    }
    if (!joined) { reps.push([r, gs, toks]); r.representative = true; r.group = r.concept; }
  }
  // a lone word is shown with the phrase that >= 90 % of its genes share
  for (const r of rows) {
    if (!r.representative || r.layer !== "word") continue;
    const cnt = new Map();
    for (const g of r.genes) for (const jp of phraseOf.get(g) || []) {
      const lab = list[jp].label;
      if (lab.split(" ").includes(r.label)) cnt.set(jp, (cnt.get(jp) || 0) + 1);
    }
    const need = Math.max(2, Math.ceil(0.9 * r.genes.length));
    let best = null;
    for (const [jp, v] of cnt) {
      if (v < need) continue;
      const a = [-v, list[jp].end - list[jp].start, list[jp].key];
      if (!best || a[0] < best[0] || (a[0] === best[0] && (a[1] < best[1] || (a[1] === best[1] && pyCmp(a[2], best[2]) < 0)))) best = a.concat(jp);
    }
    if (best) r.display = list[best[3]].label;
  }
  const members = new Map();
  for (const r of cand) { let a = members.get(r.group); if (!a) members.set(r.group, a = []); a.push(r.label); }
  const glabel = new Map(rows.map(r => [r.concept, r.label]));
  for (const r of rows) {
    r.merged = r.representative ? (members.get(r.concept) || []).slice(1).join(", ") : "";
    r.group = r.group ? glabel.get(r.group) : "";
    for (const f of ["label", "display", "group", "merged"]) r[f] = r[f].replaceAll(NBH, "-");
    r.geneIds = r.genes.map(i => ann.genes.ids[i]);
    r.symbols = r.genes.map(i => ann.genes.symbols[i] || ann.genes.ids[i]);
  }
  return { rows, study, missing, notInBackground, notAnnotated, N, n, m, params: { ...o, background_genes: bgIds.size } };
}

// ------------------------------------------------------------------------------------------- export
export function toTSV(res) {
  const cols = ["concept", "layer", "label", "k", "n", "K", "N", "fold", "p", "fdr", "log10_p", "log10_fdr", "tested",
    "status", "representative", "group", "display", "merged", "genes", "symbols"];
  const val = (r, c) => c === "genes" ? r.geneIds.join(", ") : c === "symbols" ? r.symbols.join(", ")
    : c === "representative" || c === "tested" ? (r[c] ? "True" : "False") : r[c];
  return [cols.join("\t"), ...res.rows.map(r => cols.map(c => val(r, c)).join("\t"))].join("\n") + "\n";
}

/** FASTA-like gene lists per concept (as in GeneCloud 2015, for VirtualPlant). */
export function toFasta(res) {
  return res.rows.filter(r => r.representative && r.status).map(r =>
    `>${r.display} | ${r.status} | ${r.k}/${r.n} genes | fold ${r.fold.toFixed(2)} | FDR ${r.fdr.toExponential(2)}\n` +
    r.geneIds.join("\n")).join("\n") + "\n";
}

export { LAYER_GROUP };
