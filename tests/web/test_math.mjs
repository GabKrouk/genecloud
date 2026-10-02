// Browser engine maths, no data needed: hypergeometric tail against scipy, BH / BY against direct definitions.
// usage: node tests/web/test_math.mjs
import fs from "node:fs";
import { hypergeomSF, log10HypergeomSF, log10Adjust, bh } from "../../web/engine.js";

let fail = 0;
const ref = JSON.parse(fs.readFileSync(new URL("./hypergeom_ref.json", import.meta.url)));
for (const [k, N, K, n, p, lp] of ref) {
  const v = hypergeomSF(k, N, K, n), lv = log10HypergeomSF(k, N, K, n);
  if (p > 1e-290 && Math.abs(v - p) > 1e-9 * Math.max(p, 1e-300)) { console.log("hypergeom", k, N, K, n, v, p); fail++; }
  if (Number.isFinite(lp) && Math.abs(lv - lp) > 1e-8 * Math.max(1, Math.abs(lp))) { console.log("log10", k, N, K, n, lv, lp); fail++; }
}
// BH / BY with extra hypotheses at P = 1
const p = [1e-5, 0.003, 0.04, 0.2, 0.5, 1e-9], m = 20;
const full = p.concat(Array(m - p.length).fill(1));
const order = [...full.keys()].sort((a, b) => full[a] - full[b]);
const ref_bh = new Array(m);
for (let r = 0; r < m; r++) {
  let best = Infinity;
  for (let j = r; j < m; j++) best = Math.min(best, full[order[j]] * m / (j + 1));
  ref_bh[order[r]] = Math.min(1, best);
}
let H = 0; for (let i = 1; i <= m; i++) H += 1 / i;
const q = bh(p, m), lq = log10Adjust(p.map(Math.log10), m), lqy = log10Adjust(p.map(Math.log10), m, "BY");
p.forEach((_, i) => {
  if (Math.abs(q[i] - ref_bh[i]) > 1e-12) { console.log("bh", i, q[i], ref_bh[i]); fail++; }
  if (Math.abs(10 ** lq[i] - ref_bh[i]) > 1e-12) { console.log("log bh", i); fail++; }
  if (Math.abs(10 ** lqy[i] - Math.min(1, ref_bh[i] * H)) > 1e-12) { console.log("by", i); fail++; }
});
console.log(fail ? `${fail} failure(s)` : `ok: ${ref.length} hypergeometric values, BH and BY`);
process.exit(fail ? 1 : 0);
