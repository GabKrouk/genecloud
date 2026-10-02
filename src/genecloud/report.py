"""Self-contained interactive HTML report: hover a word for its statistics, click it to see the genes."""
from __future__ import annotations

import html
import json
from pathlib import Path

from .cloud import layout

_CSS = """
:root{--bg:#fff;--fg:#111827;--sub:#6b7280;--rule:#e5e7eb;--hi:#fef3c7;--card:#f9fafb}
@media (prefers-color-scheme:dark){:root{--bg:#0b0f14;--fg:#f3f4f6;--sub:#9ca3af;--rule:#1f2937;--hi:#3b2f0b;--card:#111827}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}
.wrap{max-width:1100px;margin:0 auto;padding:24px 16px 48px}h1{font-size:22px;margin:0 0 4px}.sub{color:var(--sub);font-size:13px}
svg text{cursor:pointer;transition:opacity .15s}svg.sel text{opacity:.25}svg.sel text.on{opacity:1}
.cloud{border:1px solid var(--rule);border-radius:10px;margin:16px 0;background:var(--bg)}
#tip{position:fixed;pointer-events:none;background:var(--card);border:1px solid var(--rule);border-radius:8px;padding:8px 10px;
 font-size:12px;max-width:340px;box-shadow:0 6px 18px rgba(0,0,0,.15);display:none;z-index:9}
table{border-collapse:collapse;width:100%;font-size:12.5px}th,td{text-align:left;padding:5px 8px;border-bottom:1px solid var(--rule);vertical-align:top}
th{color:var(--sub);font-weight:600;cursor:pointer;user-select:none}tr.on td{background:var(--hi)}td.num{text-align:right;font-variant-numeric:tabular-nums}
.scroll{overflow-x:auto}h2{font-size:16px;margin:28px 0 8px}.pill{display:inline-block;border-radius:4px;padding:0 5px;font-size:11px;background:var(--card);border:1px solid var(--rule)}
#detail{background:var(--card);border:1px solid var(--rule);border-radius:10px;padding:10px 14px;min-height:44px}
"""

_JS = """
const C=__CONCEPTS__,G=__GENES__;const tip=document.getElementById('tip'),svg=document.querySelector('svg');
const f=(x,d=2)=>x<0.001?x.toExponential(1):x.toFixed(d);
svg.querySelectorAll('text').forEach(t=>{const c=C[t.dataset.i];
 t.addEventListener('mousemove',e=>{tip.style.display='block';tip.style.left=(e.clientX+14)+'px';tip.style.top=(e.clientY+14)+'px';
  tip.innerHTML=`<b>${c.display}</b> <span class=pill>${c.layer}</span><br>${c.k} of ${c.n} genes (background ${c.K} of ${c.N})<br>`+
  `fold ×${f(c.fold,1)} · P ${f(c.p)} · FDR ${f(c.fdr)}`+(c.merged?`<br><span style="color:var(--sub)">also: ${c.merged}</span>`:'')});
 t.addEventListener('mouseleave',()=>tip.style.display='none');
 t.addEventListener('click',()=>select(+t.dataset.i));});
function select(i){const c=C[i];const on=svg.classList.contains('sel')&&svg.dataset.cur==i;
 svg.classList.toggle('sel',!on);svg.dataset.cur=on?'':i;svg.querySelectorAll('text').forEach(t=>t.classList.toggle('on',!on&&+t.dataset.i===i));
 const gs=new Set(on?[]:c.genes.split(', '));document.querySelectorAll('#genes tbody tr').forEach(r=>r.classList.toggle('on',gs.has(r.dataset.g)));
 document.getElementById('detail').innerHTML=on?'<span style="color:var(--sub)">Click a word to see the genes that carry it.</span>':
  `<b>${c.display}</b> — ${c.k} genes: ${c.symbols}`+(c.merged?`<br><span style="color:var(--sub)">merged concepts: ${c.merged}</span>`:'');}
document.querySelectorAll('table').forEach(tb=>tb.querySelectorAll('th').forEach((th,j)=>th.addEventListener('click',()=>{
 const rows=[...tb.tBodies[0].rows];const num=th.classList.contains('n');const dir=th.dataset.d=th.dataset.d==='1'?'-1':'1';
 rows.sort((a,b)=>{const x=a.cells[j].dataset.v??a.cells[j].textContent,y=b.cells[j].dataset.v??b.cells[j].textContent;
  return (num?(+x-+y):x.localeCompare(y))*dir});rows.forEach(r=>tb.tBodies[0].appendChild(r));})));
"""


def write_html(result, path, gc=None, title: str = "GeneCloud", width: int = 1000, height: int = 560, **kw):
    words, norm, th = layout(result, width, height, **kw)
    concepts = [p.row for p in words]
    for c in concepts:
        c.pop("concept", None)
    svg = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" aria-label="word cloud">']
    for i, p in enumerate(words):
        style = "font-style:italic;" if p.style == "italic" else ""
        svg.append(f'<text data-i="{i}" x="{p.x:.1f}" y="{height - p.y:.1f}" text-anchor="middle" dominant-baseline="central" '
                   f'font-size="{p.size:.1f}" font-weight="{"700" if p.weight == "bold" else "400"}" fill="{p.color}" '
                   f'style="{style}font-family:DejaVu Sans,Verdana,sans-serif">{html.escape(p.text)}</text>')
    svg.append("</svg>")
    # concept table: significant + trend representatives with their merged synonyms
    t = result.table[result.table.representative]
    crow = "".join(
        f"<tr><td>{html.escape(r.display)}</td><td><span class=pill>{r.layer}</span></td><td>{r.status}</td>"
        f"<td class=num data-v={r.k}>{r.k}/{r.n}</td><td class=num data-v={r.K}>{r.K:,}</td><td class=num data-v={r.fold}>{r.fold:.1f}</td>"
        f"<td class=num data-v={r.p}>{r.p:.2g}</td><td class=num data-v={r.fdr}>{r.fdr:.2g}</td><td>{html.escape(r.symbols)}</td>"
        f"<td style='color:var(--sub)'>{html.escape(r.merged)}</td></tr>" for r in t.itertuples())
    grow = ""
    if gc is not None:
        for g in result.study:
            a = gc.genes.get(g)
            desc = ""
            if a:
                desc = (a.tair_description[0] if getattr(a, "tair_description", None) else
                        a.protein_names[0] if a.protein_names else a.description) or ""
            grow += f"<tr data-g='{g}'><td>{g}</td><td>{html.escape(gc.symbol(g))}</td><td>{html.escape(desc)}</td></tr>"
    n, N = len(result.study), result.background_size
    miss = f" · not annotated / not in background: {', '.join(result.missing)}" if result.missing else ""
    page = f"""<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>{_CSS}</style><div class=wrap>
<h1>{html.escape(title)}</h1><div class=sub>{n} genes vs {N:,} background genes{html.escape(miss)}<br>
Word size = −log10 FDR (hypergeometric test, Benjamini–Hochberg) · colour = fold enrichment · grey = trend (P ≤ {result.params.get('trend_p')}, FDR &gt; {result.params['fdr']}) · italic = UniProt keyword.
Hover a word for its statistics, click it to highlight its genes.</div>
<div class=cloud>{''.join(svg)}</div><div id=detail><span style="color:var(--sub)">Click a word to see the genes that carry it.</span></div>
<h2>Concepts</h2><div class=scroll><table><thead><tr><th>concept</th><th>layer</th><th>status</th><th class=n>genes</th><th class=n>background</th>
<th class=n>fold</th><th class=n>P</th><th class=n>FDR</th><th>genes</th><th>merged with</th></tr></thead><tbody>{crow}</tbody></table></div>
<h2>Genes</h2><div class=scroll><table id=genes><thead><tr><th>id</th><th>symbol</th><th>name</th></tr></thead><tbody>{grow}</tbody></table></div>
<p class=sub>GeneCloud 2 · annotation: {html.escape(json.dumps(result.params.get('annotation', {}).get('sources', {})))} built {result.params.get('annotation', {}).get('built', '')}</p>
<div id=tip></div><script>{_JS.replace('__CONCEPTS__', json.dumps(concepts, default=str)).replace('__GENES__', '{}')}</script></div></html>"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(page, encoding="utf-8")
    return path
