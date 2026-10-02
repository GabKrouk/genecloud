// Compare the JavaScript engine (web/engine.js) with the Python engine on the paper examples.
// usage: node tests/web/test_engine.mjs web/data REFDIR      (REFDIR: Python TSVs from tests/web/make_reference.py)
import fs from "node:fs";
import zlib from "node:zlib";
import path from "node:path";
import { Annotation, decodeIndex, run, hypergeomSF, log10HypergeomSF } from "../../web/engine.js";

const [dataDir, refDir] = process.argv.slice(2);
const gz = f => zlib.gunzipSync(fs.readFileSync(path.join(dataDir, f)));
const meta = JSON.parse(fs.readFileSync(path.join(dataDir, "meta.json")));
const genes = JSON.parse(gz("genes.json.gz"));
const blocks = {};
for (const name of Object.keys(meta.blocks)) {
  const c = JSON.parse(gz(`concepts_${name}.json.gz`));
  Object.assign(c, decodeIndex(new Uint8Array(gz(`index_${name}.bin.gz`)), c.labels.length));
  blocks[name] = c;
}
const ann = new Annotation(meta, genes, blocks);
const cases = JSON.parse(fs.readFileSync(path.join(refDir, "cases.json")));
let fail = 0;
const close = (a, b, tol) => Math.abs(a - b) <= tol * Math.max(Math.abs(a), Math.abs(b), 1e-300);
for (const cs of cases) {
  const ids = f => ann.parseIds(fs.readFileSync(f, "utf8").split("\n").filter(l => !l.startsWith("#")).join("\n"));
  const res = run(ann, ids(cs.genes), ids(cs.background), cs.opts);
  const ref = fs.readFileSync(path.join(refDir, cs.name + ".tsv"), "utf8").trim().split("\n").map(l => l.split("\t"));
  const h = ref[0], R = ref.slice(1).map(r => Object.fromEntries(h.map((c, i) => [c, r[i]])));
  const errs = [];
  if (R.length !== res.rows.length) errs.push(`rows ${res.rows.length} vs ${R.length}`);
  const byKey = new Map(res.rows.map(r => [r.concept, r]));
  let order = 0;
  R.forEach((p, i) => {
    const r = byKey.get(p.concept);
    if (!r) { errs.push(`missing ${p.concept}`); return; }
    if (+p.k !== r.k || +p.K !== r.K || +p.N !== r.N || +p.n !== r.n) errs.push(`${p.concept} counts`);
    if (Math.abs(+p.log10_p - r.log10_p) > 1e-8 * Math.max(1, Math.abs(r.log10_p))) errs.push(`${p.concept} log10 p ${r.log10_p} vs ${p.log10_p}`);
    if (Math.abs(+p.log10_fdr - r.log10_fdr) > 1e-8 * Math.max(1, Math.abs(r.log10_fdr))) errs.push(`${p.concept} log10 fdr ${r.log10_fdr} vs ${p.log10_fdr}`);
    if (p.status !== r.status) errs.push(`${p.concept} status ${r.status} vs ${p.status}`);
    if ((p.representative === "True") !== r.representative) errs.push(`${p.concept} representative`);
    if (p.display !== r.display) errs.push(`${p.concept} display '${r.display}' vs '${p.display}'`);
    if ((p.merged || "") !== r.merged) errs.push(`${p.concept} merged '${r.merged}' vs '${p.merged}'`);
    if (res.rows[i] && res.rows[i].concept !== p.concept) order++;
  });
  if (order) errs.push(`${order} rows in a different order`);
  const sig = res.rows.filter(r => r.status === "significant" && r.representative).length;
  console.log(`${errs.length ? "FAIL" : "ok  "} ${cs.name}: ${res.rows.length} rows, ${sig} significant representatives, m=${res.m}`);
  errs.slice(0, 12).forEach(e => console.log("     ", e));
  fail += errs.length > 0;
}
// hypergeometric against reference values computed by scipy
let nh = 0;
for (const [k, N, K, n, p, lp] of JSON.parse(fs.readFileSync(path.join(refDir, "hypergeom.json")))) {
  const v = hypergeomSF(k, N, K, n), lv = log10HypergeomSF(k, N, K, n);
  nh++;
  if (!close(v, p, 1e-9) && p > 1e-290) { console.log("FAIL hypergeom", k, N, K, n, v, p); fail++; }
  if (Number.isFinite(lp) && Math.abs(lv - lp) > 1e-8 * Math.max(1, Math.abs(lp))) { console.log("FAIL log10 hypergeom", k, N, K, n, lv, lp); fail++; }
}
console.log(`hypergeometric: ${nh} random cases checked against scipy`);
console.log(fail ? `${fail} failure(s)` : "all JavaScript results match Python");
process.exit(fail ? 1 : 0);
