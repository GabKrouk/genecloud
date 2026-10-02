// Visual QA of the web cloud: no overlapping words, nothing outside the picture (desktop and phone widths).
// usage: (cd web && python -m http.server 8765) & node tests/web/check_layout.mjs [http://localhost:8765/index.html]
import { chromium } from "playwright";
const URL_ = process.argv[2] || "http://localhost:8765/index.html";
const b = await chromium.launch(process.env.CHROMIUM ? { executablePath: process.env.CHROMIUM } : {});
let bad = 0;
for (const width of [1280, 390]) {
  const p = await b.newPage({ viewport: { width, height: 900 } });
  const errs = []; p.on("pageerror", e => errs.push(e.message));
  await p.goto(URL_);
  await p.waitForFunction(() => document.querySelector("#status").textContent === "", null, { timeout: 60000 });
  for (const ex of ["hormone_IAA_up", "PHR1_induced", "nitrate_responsive"]) {
    await p.click(`[data-ex=${ex}]`);
    await p.waitForFunction(() => /done/.test(document.querySelector("#status").textContent), null, { timeout: 60000 });
    for (const mode of ["fdr", "fold"]) {
      await p.check(mode === "fdr" ? "#sizeFdr" : "#sizeFold");
      await p.waitForTimeout(300);
      const r = await p.evaluate(() => {
        const svg = document.querySelector("#cloud svg"); const vb = svg.viewBox.baseVal;
        const bs = [...svg.querySelectorAll("text[data-i]")].map(t => t.getBBox());
        let ov = 0, out = 0;
        for (let i = 0; i < bs.length; i++) {
          const a = bs[i];
          if (a.x < 0 || a.y < 0 || a.x + a.width > vb.width || a.y + a.height > vb.height) out++;
          for (let j = i + 1; j < bs.length; j++) { const c = bs[j];
            const ix = Math.min(a.x + a.width, c.x + c.width) - Math.max(a.x, c.x), iy = Math.min(a.y + a.height, c.y + c.height) - Math.max(a.y, c.y);
            if (ix > 0.5 && iy > 0.5) ov++; }
        }
        return { n: bs.length, ov, out, w: vb.width };
      });
      const ok = !r.ov && !r.out; bad += !ok;
      console.log(`${ok ? "ok  " : "FAIL"} width=${width} ${ex} ${mode} words=${r.n} overlaps=${r.ov} outside=${r.out} viewBox=${r.w}`);
    }
    await p.check("#sizeFdr");
  }
  if (errs.length) { console.log("page errors:", errs); bad++; }
  await p.close();
}
await b.close(); process.exit(bad ? 1 : 0);
