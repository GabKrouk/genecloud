// Bridge with Régine (https://gabkrouk.github.io/Regine/): both are static pages of the same site.
//   GeneCloud -> Régine  the genes carrying a term (typically a GO term) -> the DAP-seq transcription factors most associated with them
//   Régine -> GeneCloud  genes bound by a TF -> their enriched words, keywords and GO terms
// Payload {genes, name, go, from}; short lists travel in the URL hash, any size by postMessage (protocol: Regine/INTEGRATION.md).
export const ORIGIN = "https://gabkrouk.github.io";
export const REGINE = ORIGIN + "/Regine/";
const AGI = /AT[1-5CM]G\d{5}/gi;
const agis = a => [...new Set((Array.isArray(a) ? a.join(" ") : String(a)).match(AGI)?.map(x => x.toUpperCase()) || [])];

export function sendToRegine({ genes, name = "", go = "" }) {
  genes = agis(genes);
  if (genes.length < 2) return false;
  const p = new URLSearchParams({ from: "genecloud", hs: "1" });
  if (name) p.set("n", name);
  if (go) p.set("go", go);
  if (genes.length <= 400) p.set("g", genes.join(","));
  const w = window.open(REGINE + "#" + p.toString(), "regine");
  if (!w) return false;
  const on = e => {
    if (e.source !== w || e.origin !== ORIGIN || e.data?.type !== "regine:ready") return;
    removeEventListener("message", on);
    w.postMessage({ type: "genecloud:genes", genes, name, go }, ORIGIN);
  };
  addEventListener("message", on);
  setTimeout(() => removeEventListener("message", on), 20000);
  return true;
}

// Lists coming from Régine: by postMessage, or in the hash (#from=regine&n=…&g=AT…,AT…)
export function incoming(onGenes) {
  addEventListener("message", e => {
    if (e.origin !== ORIGIN || e.data?.type !== "regine:genes") return;
    const genes = agis(e.data.genes || []);
    if (genes.length) onGenes({ genes, name: String(e.data.name || "") });
  });
  const p = new URLSearchParams(location.hash.replace(/^#/, ""));
  if (p.get("from") === "regine") {
    const genes = agis(p.get("g") || "");
    if (genes.length) onGenes({ genes, name: p.get("n") || "" });
    return { waiting: !genes.length && p.get("hs") === "1", name: p.get("n") || "" };
  }
  return null;
}
export const announceReady = () => { if (window.opener) window.opener.postMessage({ type: "genecloud:ready" }, ORIGIN); };
