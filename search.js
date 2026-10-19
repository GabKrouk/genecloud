// GeneCloud — search the annotation with a boolean query: which genes mention "nitrate", "xylem" …?
//
//   nitrate                      genes whose annotation contains a word starting with "nitrate" (nitrate, nitrates)
//   "root hair"                  a phrase (space, hyphen or slash between the words)
//   nitrate AND xylem            both (AND is implicit: "nitrate xylem" is the same)
//   nitrate OR auxin             either
//   nitrate NOT ammonium         the first without the second (also: nitrate -ammonium)
//   (nitrate OR nitrogen) AND (xylem OR phloem)
//   transport*                   explicit prefix; with "whole words", only words written exactly so match, except with *
//
// French operators work too in capitals (ET, OU, SAUF / NON), as do &, | and !. Case and accents are ignored.
// Precedence: NOT > AND > OR, as in PubMed / Google.

export const norm = s => String(s ?? "").normalize("NFD").replace(/[̀-ͯ]/g, "").toLowerCase();

const OPS = { and: "AND", et: "AND", "&": "AND", "&&": "AND", "+": "AND",
  or: "OR", ou: "OR", "|": "OR", "||": "OR",
  not: "NOT", non: "NOT", sauf: "NOT", "!": "NOT" };

export function tokenize(q) {
  const out = [];
  const re = /\s*(?:(\()|(\))|"([^"]*)"?|(&&|\|\||[&|!+])|(-)(?=[^\s-])|([^\s()"&|!]+))/gy;
  let m;
  q = String(q ?? "").replace(/[“”«»]/g, '"');
  while (re.lastIndex < q.length && (m = re.exec(q))) {
    if (m[1]) out.push({ t: "(" });
    else if (m[2]) out.push({ t: ")" });
    else if (m[3] !== undefined) { if (/[\p{L}\p{N}]/u.test(m[3])) out.push({ t: "term", v: m[3].trim(), phrase: true }); }
    else if (m[4]) out.push({ t: OPS[m[4]] });
    else if (m[5]) out.push({ t: "NOT" });
    else if (m[6]) {
      const w = m[6], op = OPS[w.toLowerCase()];
      // and / or / not in any case; the French ET, OU, SAUF, NON only in capitals ("non-specific" stays searchable)
      if (op && (/^(and|or|not)$/i.test(w) || w === w.toUpperCase())) out.push({ t: op });
      else if (/[\p{L}\p{N}]/u.test(w)) out.push({ t: "term", v: w });
    }
  }
  return out;
}

/** query -> tree {op:"AND"|"OR", a, b} | {op:"NOT", a} | {op:"TERM", v, phrase}; throws on a syntax error */
export function parseQuery(q) {
  const tk = tokenize(q);
  let i = 0;
  const peek = () => tk[i]?.t;
  const or = () => {
    let a = and();
    while (peek() === "OR") { i++; a = { op: "OR", a, b: and() }; }
    return a;
  };
  const and = () => {
    let a = unary();
    for (;;) {
      if (peek() === "AND") { i++; a = { op: "AND", a, b: unary() }; }
      else if (peek() === "term" || peek() === "(" || peek() === "NOT") a = { op: "AND", a, b: unary() };
      else return a;
    }
  };
  const unary = () => {
    if (peek() === "NOT") { i++; return { op: "NOT", a: unary() }; }
    if (peek() === "(") {
      i++;
      const a = or();
      if (peek() !== ")") throw new Error("a parenthesis is not closed");
      i++;
      return a;
    }
    if (peek() === "term") { const t = tk[i++]; return { op: "TERM", v: t.v, phrase: !!t.phrase }; }
    throw new Error(i < tk.length ? `unexpected ${tk[i].t === ")" ? "“)”" : tk[i].t}` : "the query ends with an operator");
  };
  if (!tk.length) return null;
  const tree = or();
  if (i < tk.length) throw new Error(tk[i].t === ")" ? "a “)” has no “(”" : `unexpected ${tk[i].t}`);
  return tree;
}

const escRe = s => s.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
const W = "a-z0-9";

/** regular expression (on normalised text) for one term */
export function termRegex(v, wholeWords = false) {
  const words = norm(v).split(/[\s\-_/]+/).filter(Boolean);
  const star = /\*$/.test(words[words.length - 1] || "");
  const body = words.map(w => w.split("*").map(escRe).join(`[${W}]*`)).join("[\\s\\-_/]+");
  const end = star || !wholeWords ? "" : `(?![${W}])`;
  return new RegExp(`(?<![${W}])${body}${end}`, "g");
}

export function terms(tree, out = [], neg = false) {
  if (!tree) return out;
  if (tree.op === "TERM") out.push({ v: tree.v, neg });
  else if (tree.op === "NOT") terms(tree.a, out, !neg);
  else { terms(tree.a, out, neg); terms(tree.b, out, neg); }
  return out;
}

/** compile the query: returns {test(text), positive: [RegExp…]}; text must be normalised with norm() */
export function compile(q, { wholeWords = false } = {}) {
  const tree = parseQuery(q);
  if (!tree) return null;
  const cache = new Map();
  const re = v => cache.get(v) || cache.set(v, termRegex(v, wholeWords)).get(v);
  const ev = (n, s) => {
    switch (n.op) {
      case "TERM": { const r = re(n.v); r.lastIndex = 0; return r.test(s); }
      case "NOT": return !ev(n.a, s);
      case "AND": return ev(n.a, s) && ev(n.b, s);
      default: return ev(n.a, s) || ev(n.b, s);
    }
  };
  const pos = [...new Set(terms(tree).filter(t => !t.neg).map(t => t.v))].map(v => termRegex(v, wholeWords));
  return { tree, test: s => ev(tree, s), positive: pos };
}

/** text -> [[start, end]…] of the positive terms, on the original text (normalisation keeps the length of plain letters) */
export function matches(text, positive) {
  const n = norm(text);
  if (n.length !== String(text).length) return [];      // rare composed characters: no highlight rather than a wrong one
  const out = [];
  for (const r of positive) {
    r.lastIndex = 0;
    for (const m of n.matchAll(r)) {
      let e = m.index + m[0].length;
      while (e < n.length && /[a-z0-9]/.test(n[e])) e++;   // highlight the whole word (nitrate|s)
      out.push([m.index, e]);
    }
  }
  out.sort((a, b) => a[0] - b[0]);
  const merged = [];
  for (const x of out) if (merged.length && x[0] <= merged[merged.length - 1][1]) merged[merged.length - 1][1] = Math.max(merged[merged.length - 1][1], x[1]); else merged.push(x);
  return merged;
}
