// GeneCloud web application: loads the annotation files, runs the engine, draws the cloud.
import { Annotation, decodeIndex, run, toTSV, toFasta } from "./engine.js";
import { sendToRegine, incoming, announceReady } from "./regine.js";

const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];
const DATA = "data/";

// ------------------------------------------------------------------------------------------ banner
(() => {
  let s = "";
  for (let r = 0; r < 8; r++) {
    let line = "";
    for (let c = 0; c < 90; c++) line += Math.random() < 0.18 ? (Math.random() < 0.5 ? "0" : "1") : " ";
    s += line + "\n";
  }
  $("#bits").textContent = s;
})();
$("#theme").addEventListener("click", () => {
  const root = document.documentElement;
  const dark = root.dataset.theme ? root.dataset.theme === "dark"
    : matchMedia("(prefers-color-scheme: dark)").matches;
  root.dataset.theme = dark ? "light" : "dark";
  try { localStorage.setItem("gc-theme-2026", root.dataset.theme); } catch (e) {}
});
try { const t = localStorage.getItem("gc-theme-2026"); if (t) document.documentElement.dataset.theme = t; } catch (e) {}

// ------------------------------------------------------------------------------------------ loading
async function fetchBytes(f) {
  const r = await fetch(DATA + f);
  if (!r.ok) throw new Error(`cannot load ${f} (${r.status})`);
  let buf = new Uint8Array(await r.arrayBuffer());
  if (buf[0] === 0x1f && buf[1] === 0x8b) {          // gzip (served as a file, not as Content-Encoding)
    const ds = new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"));
    buf = new Uint8Array(await new Response(ds).arrayBuffer());
  }
  return buf;
}
const fetchJSON = async f => JSON.parse(new TextDecoder().decode(await fetchBytes(f)));

let ANN = null, META = null, TEXTS = null;
const blockPromises = {};
function loadBlock(name) {
  return blockPromises[name] ??= (async () => {
    const [c, idx] = await Promise.all([fetchJSON(`concepts_${name}.json.gz`), fetchBytes(`index_${name}.bin.gz`)]);
    Object.assign(c, decodeIndex(idx, c.labels.length));
    ANN.blocks[name] = c;
    ANN._cache?.clear();
  })();
}
const status = (html) => { $("#status").innerHTML = html; };
const ready = (async () => {
  status('<span class="loader"></span>loading the annotation…');
  META = await (await fetch(DATA + "meta.json")).json();
  const genes = await fetchJSON("genes.json.gz");
  ANN = new Annotation(META, genes, {});
  await Promise.all([loadBlock("wp_all"), loadBlock("kw")]);
  for (const el of $$("[data-n]")) {
    const n = META.background_sizes[el.dataset.n];
    if (n) el.textContent = `(${n.toLocaleString("en")} annotated genes)`;
  }
  const a = META.annotation;
  const rel = (a.tair_release || "").replace("TAIR_Data_", "");
  $("#annot").innerHTML = `Annotation built ${a.built}: ${META.n_genes.toLocaleString("en")} genes · TAIR public release `
    + `${rel ? rel.slice(0, 4) + "-" + rel.slice(4, 6) + "-" + rel.slice(6) : "—"} (Araport11 descriptions, curator summaries, phenotypes, symbols) · `
    + `UniProtKB reference proteome · NCBI Gene · Gene Ontology.`;
  status("");
  return ANN;
})().catch(e => { status(`<span class="err">${e.message}</span>`); throw e; });

// ------------------------------------------------------------------------------------------- form
const countIds = () => {
  const n = ANN ? ANN.parseIds($("#genes").value).length : ($("#genes").value.match(/AT[1-5CM]G\d{5}/gi) || []).length;
  $("#gcount").textContent = `${n} identifier${n === 1 ? "" : "s"}`;
};
$("#genes").addEventListener("input", countIds);
$$("input[name=bg]").forEach(r => r.addEventListener("change", () => {
  $("#bgcustom").style.display = r.value === "custom" && r.checked ? "block" : "none";
  $("#bgcount").textContent = "";
}));
$("#bgcustom").addEventListener("input", () => {
  const n = ANN ? ANN.parseIds($("#bgcustom").value).length : 0;
  $("#bgcount").textContent = `${n} identifiers`;
});
$$(".examples button").forEach(b => b.addEventListener("click", async () => {
  const t = await (await fetch(`examples/${b.dataset.ex}.txt`)).text();
  $("#genes").value = t.split("\n").filter(l => !l.startsWith("#")).join("\n").trim();
  const ath1 = $("input[name=bg][value=ATH1]"); ath1.checked = true; ath1.dispatchEvent(new Event("change"));
  countIds();
  $("#form").requestSubmit();
}));
$("#reset").addEventListener("click", () => {
  setTimeout(() => { countIds(); $("#bgcustom").style.display = "none"; $("#results").style.display = "none"; status(""); history.replaceState(null, "", location.pathname); });
});

function options() {
  const v = n => $(`input[name=${n}]:checked`)?.value;
  return { preset: v("preset"), layers: $$("input[name=layer]:checked").map(x => x.value),
    minGenes: +v("min"), fdr: +v("fdr"), adjust: v("adjust") };
}

function background() {
  const choice = $("input[name=bg]:checked")?.value;
  if (!choice) return null;
  if (choice === "custom") return ANN.parseIds($("#bgcustom").value);
  const bit = META.backgrounds[choice];
  return ANN.genes.ids.filter((g, i) => ANN.genes.flags[i] & bit);
}

let RES = null, LAST = null;
$("#form").addEventListener("submit", async e => {
  e.preventDefault();
  $("#submit").disabled = true;
  try {
    await ready;
    const genes = ANN.parseIds($("#genes").value);
    if (!genes.length) throw new Error("Paste at least one AGI identifier (e.g. AT1G01010).");
    const bg = background();
    if (!bg) throw new Error("Choose a background: GeneCloud needs to know which genes could have been in your list.");
    if (!bg.length) throw new Error("The background list contains no AGI identifier.");
    const o = options();
    if (!o.layers.length) throw new Error("Select at least one kind of concept.");
    status('<span class="loader"></span>computing…');
    const need = [];
    if (o.layers.includes("word") || o.layers.includes("phrase")) need.push("wp_" + o.preset);
    if (o.layers.includes("keyword")) need.push("kw");
    if (o.layers.includes("go")) need.push("go");
    await Promise.all(need.map(loadBlock));
    await new Promise(r => setTimeout(r, 20));
    const t0 = performance.now();
    RES = run(ANN, genes, bg, o);
    LAST = { genes, o, bgChoice: $("input[name=bg]:checked").value };
    show(RES);
    status(`done in ${((performance.now() - t0) / 1000).toFixed(1)} s`);
    $("#results").scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (err) {
    status(`<span class="err">${err.message}</span>`);
  } finally {
    $("#submit").disabled = false;
  }
});

// ------------------------------------------------------------------------------------------ cloud
const RAMP = { light: ["#8fc9c4", "#2f8f9d", "#2a5ea8", "#5b3a9e", "#b0216f"],
               dark: ["#6b7280", "#38bdf8", "#facc15", "#f97316", "#ef4444"] };
const hex = h => [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
function ramp(stops, t) {
  t = Math.max(0, Math.min(1, t));
  const x = t * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(x)), f = x - i;
  const a = hex(stops[i]), b = hex(stops[i + 1]);
  return "#" + a.map((v, j) => Math.round(v + (b[j] - v) * f).toString(16).padStart(2, "0")).join("");
}
const FONT = '"Helvetica Neue",Helvetica,Arial,sans-serif';
const fmt = (x, d = 2) => x === 0 ? "0" : x < 1e-3 ? x.toExponential(d === 2 ? 1 : d) : x.toPrecision(d === 2 ? 2 : 3);
const sup = n => String(n).replace(/[0-9-]/g, c => "⁰¹²³⁴⁵⁶⁷⁸⁹"["0123456789".indexOf(c)] ?? "⁻");
const esc = s => String(s).replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

/** deterministic spiral layout, most significant first; shrinks until every significant term fits */
function layout(res, W, H, opt) {
  const sig = res.rows.filter(r => r.representative && r.status === "significant").slice(0, opt.top);
  const tr = opt.showTrend ? res.rows.filter(r => r.representative && r.status === "trend").slice(0, opt.maxTrend) : [];
  const ctx = document.createElement("canvas").getContext("2d");
  const lf = sig.map(r => Math.log2(r.fold));
  const lo = Math.max(0, Math.floor(Math.min(...lf, 1))), hi = Math.max(lo + 1, Math.ceil(Math.max(...lf, 1)));
  const byFold = opt.sizeBy === "fold";
  const s = sig.map(r => byFold ? Math.log2(r.fold) : -r.log10_fdr);
  const smin = byFold ? Math.min(...s) : Math.min(...s, -Math.log10(res.params.fdr)), smax = Math.max(...s, smin);
  const label = r => opt.counts ? `${r.display}(${r.k}|${Math.round(r.fold)})` : r.display;
  let scale = 1, placed, sizeOf;
  for (let it = 0; it < 12; it++) {
    const top = Math.max(opt.minFont * 1.25, opt.maxFont * scale);
    sizeOf = v => opt.minFont + (top - opt.minFont) * (smax > smin ? (v - smin) / (smax - smin) : 1);
    const items = sig.map((r, i) => ({ r, text: label(r), fs: sizeOf(s[i]), bold: true,
      color: ramp(RAMP[opt.theme], (lf[i] - lo) / (hi - lo)) }))
      .concat(tr.map(r => ({ r, text: label(r), fs: opt.minFont * 0.95, bold: false, color: opt.theme === "dark" ? "#4b5563" : "#c3c8ce" })));
    const boxes = [];
    placed = [];
    const ax = Math.sqrt(W / H);
    for (const it2 of items) {
      ctx.font = `${it2.r.layer === "keyword" ? "italic " : ""}${it2.bold ? 700 : 400} ${it2.fs}px ${FONT}`;
      const mt = ctx.measureText(it2.text);
      const asc = mt.fontBoundingBoxAscent ?? it2.fs * 0.95, desc = mt.fontBoundingBoxDescent ?? it2.fs * 0.25;
      const w = mt.width + 4, h = asc + desc + 2;
      let pos = null;
      for (let th = 0; th < 420; th += 0.32 / (1 + th / 60)) {
        const rr = 2.1 * th;
        const x = W / 2 + rr * Math.cos(th) * ax - w / 2, y = H / 2 + rr * Math.sin(th) / ax - h / 2;
        if (x < 4 || y < 4 || x + w > W - 4 || y + h > H - 4) continue;
        if (boxes.every(([a, b, c, d]) => x + w <= a || a + c <= x || y + h <= b || b + d <= y)) { pos = [x, y, w, h]; break; }
      }
      if (pos) { boxes.push(pos); placed.push({ ...it2, x: pos[0] + w / 2, y: pos[1] + 1 + asc }); }   // y = baseline
      else if (it2.bold) { placed = null; break; }
    }
    if (placed) break;
    scale *= 0.86;
  }
  let key = [];
  if (sig.length && byFold) {
    const marks = [];
    for (let v = Math.ceil(smin); v <= Math.floor(smax); v++) marks.push(v);
    if (!marks.length) marks.push(smin);
    key = marks.filter((v, i) => marks.length <= 4 || i % Math.ceil(marks.length / 4) === 0 || i === marks.length - 1)
      .map(v => ({ fs: sizeOf(v), lab: `×${+(2 ** v).toPrecision(2)}` }));
  } else if (sig.length) {
    const marks = [smin];
    const topMark = Math.floor(smax);
    if (topMark - smin > 4) marks.push(Math.round((smin + topMark) / 2));
    if (topMark - smin > 1.5) marks.push(topMark);
    key = marks.map(v => ({ fs: sizeOf(v), lab: Math.abs(v + Math.log10(res.params.fdr)) < 1e-9 ? String(res.params.fdr) : `10${sup(-v)}` }));
  }
  return { placed: placed || [], lo, hi, key, nSig: sig.length };
}

function cloudSVG(res, opt) {
  const narrow = ($("#cloud").clientWidth || 1000) < 640;
  const W = narrow ? 640 : 1000, H = narrow ? 720 : 560;
  const L = layout(res, W, H, opt);
  const bg = opt.theme === "dark" ? "#0b0f14" : "#ffffff";
  const parts = [`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="GeneCloud word cloud">`,
    `<rect width="100%" height="100%" fill="${bg}"/>`];
  L.placed.forEach((p, i) => {
    parts.push(`<text data-i="${res.rows.indexOf(p.r)}" x="${p.x.toFixed(1)}" y="${p.y.toFixed(1)}" text-anchor="middle"`
      + ` font-family='${FONT}' font-size="${p.fs.toFixed(1)}" font-weight="${p.bold ? 700 : 400}"${p.r.layer === "keyword" ? ' font-style="italic"' : ""} fill="${p.color}">${esc(p.text)}</text>`);
  });
  if (!L.nSig) parts.push(`<text x="${W / 2}" y="${H - 24}" text-anchor="middle" font-family='${FONT}' font-size="15" fill="#8a929a">no term passes FDR ${res.params.fdr}</text>`);
  parts.push("</svg>");
  return { svg: parts.join(""), L };
}

// ------------------------------------------------------------------------------------------ display
function show(res) {
  $("#results").style.display = "block";
  const rep = res.rows.filter(r => r.representative && r.status);
  const nsig = rep.filter(r => r.status === "significant").length;
  $("#summary").innerHTML = `<b>${res.n}</b> of your ${res.n + res.missing.length} genes analysed against <b>${res.N.toLocaleString("en")}</b>`
    + ` annotated background genes · <b>${nsig}</b> significant term${nsig === 1 ? "" : "s"} (FDR ≤ ${res.params.fdr}, `
    + `${res.params.adjust === "BY" ? "Benjamini–Yekutieli" : "Benjamini–Hochberg"} over ${res.m.toLocaleString("en")} testable terms)`
    + (res.notInBackground.length ? ` · ${res.notInBackground.length} not in the background` : "")
    + (res.notAnnotated.length ? ` · ${res.notAnnotated.length} without annotation` : "");
  drawCloud();
  // table (as in Fig. 1D of the 2015 paper)
  const tb = $("#terms tbody");
  tb.innerHTML = rep.map(r => `<tr class="${r.status}" data-i="${res.rows.indexOf(r)}"><td><a href="#" class="term">${esc(r.display)}</a> <a href="#" class="rg" title="Send the ${r.k} genes of this term to Régine: which transcription factors bind them?">&rarr; Régine</a></td>`
    + `<td><span class="pill">${esc(r.layer)}</span></td><td class="num">${r.k}</td><td class="num">${r.K.toLocaleString("en")}</td>`
    + `<td class="num">${r.fold.toFixed(1)}</td><td class="num">${fmt(r.p, 3)}</td><td class="num">${fmt(r.fdr, 3)}</td>`
    + `<td class="note">${esc(r.merged)}</td></tr>`).join("");
  tb.querySelectorAll("a.rg").forEach(a => a.addEventListener("click", ev => { ev.preventDefault(); toRegine(RES.rows[+a.closest("tr").dataset.i]); }));
  tb.querySelectorAll("a.term").forEach(a => a.addEventListener("click", ev => { ev.preventDefault(); select(+a.closest("tr").dataset.i, true); }));
  // genes
  $("#glist tbody").innerHTML = res.study.map(g => {
    const i = ANN.index.get(g);
    return `<tr data-g="${g}"><td>${g}</td><td>${esc(ANN.genes.symbols[i] || "")}</td><td>${esc(ANN.genes.desc[i] || "")}</td></tr>`;
  }).join("");
  $("#missing").textContent = [res.notInBackground.length ? `Not in the background: ${res.notInBackground.join(", ")}.` : "",
    res.notAnnotated.length ? `In the background but without annotation for the chosen concepts: ${res.notAnnotated.join(", ")}.` : ""].join(" ");
  $("#locus").innerHTML = '<p class="note">Click a word in the cloud or in the table.</p>';
  tab("ptable");
}

function cloudOpts() {
  return { top: 70, maxTrend: 25, showTrend: $("#showTrend").checked, counts: $("#showCounts").checked,
    theme: $("#darkCloud").checked ? "dark" : "light", maxFont: 58, minFont: 11,
    sizeBy: $("input[name=sizeby]:checked")?.value || "fdr" };
}
let CLOUD = null;
function drawCloud() {
  if (!RES) return;
  const opt = cloudOpts();
  CLOUD = cloudSVG(RES, opt);
  $("#cloud").innerHTML = CLOUD.svg;
  $("#cloud").style.background = opt.theme === "dark" ? "#0b0f14" : "#fff";
  const { lo, hi, key } = CLOUD.L;
  const stops = RAMP[opt.theme];
  $("#ramp").style.background = `linear-gradient(90deg,${stops.join(",")})`;
  const ticks = [];
  for (let t = lo; t <= hi; t++) ticks.push(`<span>×${2 ** t}</span>`);
  $("#rampticks").innerHTML = ticks.join("");
  $("#sizekey").innerHTML = key.length ? (opt.sizeBy === "fold" ? "word size (log2 fold enrichment): " : "word size (−log10 FDR): ") + key.map(k => `<span style="font-size:${Math.min(k.fs, 30)}px">${k.lab}</span>`).join("")
    + (opt.showTrend ? ` · <span style="font-weight:400">grey = nominal trend (uncorrected P ≤ ${RES.params.trendP})</span>` : "") + ' · <i>italic</i> = UniProt keyword' : "";
  const svg = $("#cloud svg");
  svg.querySelectorAll("text[data-i]").forEach(t => {
    const r = RES.rows[+t.dataset.i];
    t.addEventListener("mousemove", e => {
      const tip = $("#tip");
      tip.style.display = "block";
      tip.style.left = Math.min(e.clientX + 14, innerWidth - 370) + "px";
      tip.style.top = (e.clientY + 14) + "px";
      tip.innerHTML = `<b>${esc(r.display)}</b> <span class="pill">${r.layer}</span><br>${r.k} of ${r.n} genes · background ${r.K} of ${r.N}`
        + `<br>fold ×${r.fold.toFixed(1)} · P ${fmt(r.p)} · FDR ${fmt(r.fdr)}`
        + (r.merged ? `<br><span class="note">also: ${esc(r.merged)}</span>` : "") + `<br><span class="note">click to see the genes</span>`;
    });
    t.addEventListener("mouseleave", () => { $("#tip").style.display = "none"; });
    t.addEventListener("click", () => select(+t.dataset.i, true));
  });
}
["#showTrend", "#showCounts", "#darkCloud", "#sizeFdr", "#sizeFold"].forEach(s => $(s).addEventListener("change", drawCloud));
let RESIZE = 0;
addEventListener("resize", () => { clearTimeout(RESIZE); RESIZE = setTimeout(drawCloud, 250); });

function tab(id) {
  $$(".tabs button").forEach(b => b.classList.toggle("on", b.dataset.p === id));
  $$(".panel").forEach(p => p.classList.toggle("on", p.id === id));
}
$$(".tabs button").forEach(b => b.addEventListener("click", () => tab(b.dataset.p)));

// genes carrying a term, with their annotation and the term highlighted (Fig. 1E of the 2015 paper)
let SELECTED = -1;
async function select(i, open) {
  const r = RES.rows[i];
  SELECTED = i;
  const svg = $("#cloud svg");
  svg.classList.add("sel");
  svg.querySelectorAll("text[data-i]").forEach(t => t.classList.toggle("on", +t.dataset.i === i));
  const gs = new Set(r.geneIds);
  $$("#glist tbody tr").forEach(tr => tr.classList.toggle("on", gs.has(tr.dataset.g)));
  $$("#terms tbody tr").forEach(tr => tr.classList.toggle("on", +tr.dataset.i === i));
  if (open) tab("plocus");
  $("#locus").innerHTML = '<p class="note"><span class="loader"></span>loading the annotation texts…</p>';
  TEXTS ??= await fetchJSON("texts.json.gz");
  const words = [r.label, ...r.merged.split(", ")].filter(Boolean).join(" ").toLowerCase()
    .split(/[\s\-]+/).filter(w => w.length > 2).map(w => w.replace(/(ies|es|s)$/, "").slice(0, Math.max(4, w.length - 2)));
  const re = words.length ? new RegExp(`\\b(${[...new Set(words)].map(w => w.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("|")})[\\w-]*`, "gi") : null;
  const hl = s => { const e = esc(s || ""); return re ? e.replace(re, m => `<mark>${m}</mark>`) : e; };
  const rows = r.geneIds.map(g => {
    const k = ANN.index.get(g), t = TEXTS[k] || [];
    return `<tr><td><a href="https://www.arabidopsis.org/locus?name=${g}" target="_blank" rel="noopener">${g}</a><br><b>${esc(ANN.genes.symbols[k] || "")}</b></td>`
      + `<td>${esc((t[5] || "").replace(/_/g, " "))}</td><td class="long">${hl(t[0])}</td><td class="long">${hl(t[1])}</td>`
      + `<td class="long">${hl(t[2])}</td><td class="long">${hl(t[3])}</td><td class="long">${hl(t[4])}</td></tr>`;
  }).join("");
  $("#locus").innerHTML = `<h3 style="margin-top:6px">Find locus containing the term ‘${esc(r.display)}’</h3>`
    + `<p class="note">${r.k} of ${r.n} genes (background ${r.K} of ${r.N}) · fold ×${r.fold.toFixed(1)} · P ${fmt(r.p)} · FDR ${fmt(r.fdr)}`
    + (r.merged ? ` · merged terms: ${esc(r.merged)}` : "") + ` · <a href="#" id="locus-rg">&rarr; Régine: transcription factors binding these ${r.k} genes</a></p><div class="scroll tall"><table><thead><tr><th>Locus</th><th>Type</th>`
    + `<th>Short description</th><th>Curator summary</th><th>Computational description</th><th>UniProt</th><th>Phenotypes</th></tr></thead>`
    + `<tbody>${rows}</tbody></table></div>`;
}

$("#locus").addEventListener("click", ev => {
  if (ev.target.id !== "locus-rg") return;
  ev.preventDefault();
  toRegine(RES.rows[SELECTED]);
});

// sortable term table
$$("#terms thead th[data-s]").forEach(th => th.addEventListener("click", () => {
  const key = th.dataset.s, num = th.classList.contains("n");
  const dir = th.dataset.d = th.dataset.d === "1" ? "-1" : "1";
  const tb = $("#terms tbody");
  [...tb.rows].sort((a, b) => {
    const x = RES.rows[+a.dataset.i][key], y = RES.rows[+b.dataset.i][key];
    return (num ? x - y : String(x).localeCompare(String(y))) * dir;
  }).forEach(tr => tb.appendChild(tr));
}));

// ------------------------------------------------------------------------------------------ export
function download(name, data, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([data], { type }));
  a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
$("#dlsvg").addEventListener("click", () => download("genecloud.svg", CLOUD.svg, "image/svg+xml"));
$("#dltsv").addEventListener("click", () => download("genecloud.tsv", toTSV(RES), "text/tab-separated-values"));
$("#dlfa").addEventListener("click", () => download("genecloud_terms.txt", toFasta(RES), "text/plain"));
$("#dlpng").addEventListener("click", () => {
  const img = new Image();
  img.onload = () => {
    const c = document.createElement("canvas"); c.width = 2000; c.height = 1120;
    c.getContext("2d").drawImage(img, 0, 0, 2000, 1120);
    c.toBlob(b => download("genecloud.png", b, "image/png"));
  };
  img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(CLOUD.svg);
});
$("#share").addEventListener("click", async () => {
  const p = new URLSearchParams();
  p.set("genes", LAST.genes.join(","));
  p.set("bg", LAST.bgChoice === "custom" ? "custom" : LAST.bgChoice);
  p.set("text", LAST.o.preset); p.set("layers", LAST.o.layers.join(",")); p.set("min", LAST.o.minGenes); p.set("fdr", LAST.o.fdr);
  p.set("adjust", LAST.o.adjust);
  const url = location.origin + location.pathname + "#" + p.toString();
  history.replaceState(null, "", url);
  try { await navigator.clipboard.writeText(url); status("link copied" + (LAST.bgChoice === "custom" ? " (your own background is not included in the link)" : "")); }
  catch (e) { status("link in the address bar"); }
});

// restore an analysis from the URL (#genes=...&bg=ATH1...)
(async () => {
  if (!location.hash.includes("genes=")) return;
  const p = new URLSearchParams(location.hash.slice(1));
  $("#genes").value = (p.get("genes") || "").split(",").join("\n");
  const set = (n, v) => { const el = $(`input[name=${n}][value="${v}"]`); if (el) el.checked = true; };
  if (p.get("bg") && p.get("bg") !== "custom") set("bg", p.get("bg"));
  if (p.get("text")) set("preset", p.get("text"));
  if (p.get("layers")) { const L = p.get("layers").split(","); $$("input[name=layer]").forEach(x => x.checked = L.includes(x.value)); }
  if (p.get("min")) set("min", p.get("min"));
  if (p.get("fdr")) set("fdr", p.get("fdr"));
  if (p.get("adjust")) set("adjust", p.get("adjust"));
  await ready; countIds();
  if ($("input[name=bg]:checked")) $("#form").requestSubmit();
})();
ready.then(countIds);

// ------------------------------------------------------------------------------------------ Régine bridge
function toRegine(r) {
  const go = /^g:GO:\d+$/.test(r.concept) ? r.concept.slice(2) : "";     // GO terms keep their identifier
  const name = r.display || r.label;
  if (!sendToRegine({ genes: r.geneIds, name: go ? name : `${name} (${r.layer})`, go })) status('<span class="err">Could not open Régine (pop-up blocked?).</span>');
}
$("#toregine").addEventListener("click", () => {
  if (!LAST || !sendToRegine({ genes: LAST.genes, name: "GeneCloud list" })) status('<span class="err">Run an analysis first.</span>');
});
async function fromRegine({ genes, name }) {     // a list sent by Régine: look at its GO terms and words
  await ready;
  $("#genes").value = genes.join("\n");
  const ath1 = $("input[name=bg][value=ATH1]"); if (ath1) { ath1.checked = true; ath1.dispatchEvent(new Event("change")); }
  $$("input[name=layer]").forEach(x => { if (x.value === "go") x.checked = true; });
  countIds();
  status(`List received from Régine${name ? ": " + esc(name) : ""} (${genes.length} genes)`);
  $("#form").requestSubmit();
}
incoming(fromRegine);
ready.then(announceReady);
