/* IonWeaver · the page. Plain JavaScript, no build step. Every number and every structure comes
   from the engine server; this file only asks, draws and filters. */
"use strict";
const API = window.IONWEAVER_API || "";
const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmt = (x, d = 0) => (x === null || x === undefined || x === "" ? "–" : Number(x).toLocaleString("en", { maximumFractionDigits: d, minimumFractionDigits: d }));
const pct = (x, d = 1) => (x === null || x === undefined ? "–" : (Number(x)).toFixed(d) + "%");
const GROUP = new Set(["Graft", "Prune", "Swap", "Cross-over"]);

const S = { boot: null, screen: "start", family: null, start: null, seed: null, gen: null, chain: null,
            sheet: null, batch: null, filters: null, sort: "order", bench: "generate" };

/* ---------------- plumbing ---------------- */
async function api(path, body) {
  const r = await fetch(API + path, body === undefined ? {} : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
  if (!r.ok) {
    let msg = r.statusText;
    try { const j = await r.json(); const d = j.detail; msg = Array.isArray(d) ? d.map((e) => `${(e.loc || []).slice(1).join(".")}: ${e.msg}`).join("; ") : (typeof d === "string" ? d : JSON.stringify(j)); } catch (_) {}
    throw new Error(msg);
  }
  return r.json();
}
async function job(path, body, onProgress) {
  let st = await api(path, body);
  while (st.state === "queued" || st.state === "running") {
    if (onProgress) onProgress(st);
    await new Promise((res) => setTimeout(res, 450));
    st = await api(`/api/job/${st.job_id}`);
  }
  if (st.state === "failed") throw new Error(st.error || "The job failed.");
  if (onProgress) onProgress(st);
  return st.result;
}
function toast(msg, ms = 2600) { const t = $("#toast"); t.textContent = msg; t.hidden = false; clearTimeout(t._t); t._t = setTimeout(() => (t.hidden = true), ms); }
function download(name, text, type = "text/plain") { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 2000); }
async function exportProducts(products, fmt, name) {
  const slim = products.map((p) => { const q = { ...p }; delete q.svg; return q; });
  const r = await fetch(API + "/api/export", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ products: slim, fmt, name }) });
  if (!r.ok) { toast("Export failed"); return; }
  download(`${name}.${fmt}`, await r.text(), fmt === "csv" ? "text/csv" : "text/plain");
  toast(`${products.length} ions saved as ${fmt.toUpperCase()}`);
}
function progressBar(el, st) {
  const bar = $(".progress", el), msg = $(".pmsg", el);
  if (bar) { bar.classList.toggle("busy", st.state === "running" && st.progress < 0.05); $("i", bar).style.width = Math.max(3, Math.round(100 * (st.progress || 0))) + "%"; }
  if (msg) msg.textContent = st.message || st.state;
}
const opchip = (op) => `<span class="opchip ${GROUP.has(op) ? "group" : ""}">${esc(op)}</span>`;
const ionType = (t) => `<span class="ion-type ${t}">${t}</span>`;
function rangeKnob(id, label, min, max, value, step, hint, unit = "") {
  return `<div class="knob"><div class="top"><label for="${id}">${label}</label><output id="${id}-out">${value}${unit}</output></div>
    <input type="range" id="${id}" min="${min}" max="${max}" step="${step}" value="${value}">${hint ? `<div class="hint">${hint}</div>` : ""}</div>`;
}
/* rules the user sets BEFORE generating: the engine keeps examining until N ions obey them.
   Body rules speak about what lies outside the protected core; bounds speak about the whole ion. */
const BOUND_KEYS = ["mw", "heavy_atoms", "hba", "hbd", "logp", "tpsa", "rotatable", "rings"];
const HALOGENS = ["F", "Cl", "Br", "I"];
function constraintsPanel(id) {
  const els = S.boot.elements || ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"];
  return `<details class="rules" id="${id}"><summary><span>Rules for the new ions</span><span class="rsum small">none · every certified ion counts</span></summary>
  <div class="rules-body">
    <div class="rgroup"><div class="rlabel">Outside the core, allow only these elements <span class="small">· <a href="#" data-q="all">all</a> · <a href="#" data-q="nohal">no halogens</a> · <a href="#" data-q="conly">carbon only</a> · <a href="#" data-q="cno">C, N, O</a></span></div>
      <div class="elems rel">${els.map((e) => `<button type="button" class="elem on" data-el="${e}">${e}</button>`).join("")}</div>
      <div class="small" style="margin-top:4px">A switched-off element is kept out of everything the edits may place: library fragments, tails and single atoms alike.</div></div>
    <div class="rgroup"><div class="rlabel">Outside the core</div>
      <div class="switches"><label class="sw"><input type="checkbox" data-r="body_multiple_bonds"> no double or triple bonds</label><label class="sw"><input type="checkbox" data-r="body_rings"> no rings</label><label class="sw"><input type="checkbox" data-r="body_aromatic_rings"> no aromatic rings</label><label class="sw"><input type="checkbox" data-r="body_heteroatoms"> no heteroatoms</label>
      <label class="sw">at most <input type="number" class="rnum" data-r="body_heavy_atoms" min="0" step="1" placeholder="–"> atoms, hydrogen not counted</label></div></div>
    <div class="rgroup"><div class="rlabel">Whole ion, between</div>
      <div class="cons-grid">${BOUND_KEYS.map((k) => { const d = S.boot.descriptors.find((x) => x.key === k); return `<div class="c" data-k="${k}"><label>${esc(d.label)}${d.unit ? " (" + d.unit + ")" : ""}</label><div class="range2"><input type="number" class="clo" placeholder="min" step="any"><input type="number" class="chi" placeholder="max" step="any"></div></div>`; }).join("")}</div></div>
    <div class="actions" style="margin-top:6px"><button class="btn sec sm cclear" type="button">Clear every rule</button><span class="small">The core itself is never touched; these rules shape only what the edits may build around it.</span></div>
  </div></details>`;
}
function readConstraints(id) {
  const root = $(`#${id}`); if (!root) return null;
  const out = {};
  const els = $$(".rel .elem", root), on = els.filter((b) => b.classList.contains("on")).map((b) => b.dataset.el);
  if (on.length < els.length) out.body_elements = on;
  $$("input[type=checkbox][data-r]", root).forEach((c) => { if (c.checked) out[c.dataset.r] = [null, 0]; });
  $$(".rnum[data-r]", root).forEach((n) => { if (n.value !== "") out[n.dataset.r] = [null, Number(n.value)]; });
  $$(".c", root).forEach((box) => { const lo = $(".clo", box).value, hi = $(".chi", box).value; if (lo !== "" || hi !== "") out[box.dataset.k] = [lo === "" ? null : Number(lo), hi === "" ? null : Number(hi)]; });
  return Object.keys(out).length ? out : null;
}
function wireConstraints(id) {
  const root = $(`#${id}`); if (!root) return;
  const sum = () => { const c = readConstraints(id); $(".rsum", root).textContent = c ? constraintsText(c) : "none · every certified ion counts"; };
  $$(".rel .elem", root).forEach((b) => b.addEventListener("click", () => { b.classList.toggle("on"); sum(); }));
  $$("[data-q]", root).forEach((a) => a.addEventListener("click", (e) => { e.preventDefault(); const q = a.dataset.q; $$(".rel .elem", root).forEach((b) => { const el = b.dataset.el; b.classList.toggle("on", q === "all" || (q === "nohal" && !HALOGENS.includes(el)) || (q === "conly" && el === "C") || (q === "cno" && ["C", "N", "O"].includes(el))); }); sum(); }));
  $$("input", root).forEach((i) => i.addEventListener("input", sum));
  $(".cclear", root).addEventListener("click", () => { $$("input", root).forEach((i) => { if (i.type === "checkbox") i.checked = false; else i.value = ""; }); $$(".rel .elem", root).forEach((b) => b.classList.add("on")); sum(); });
}
function constraintsText(c) {
  if (!c) return "";
  const parts = [];
  for (const [k, v] of Object.entries(c)) {
    if (k === "body_elements") { parts.push(`only ${v.join(", ")} outside the core`); continue; }
    const d = S.boot.descriptors.find((x) => x.key === k); const [lo, hi] = v;
    if (d && d.body && lo === null && hi === 0) { parts.push(d.label.toLowerCase().replace(/^/, "no ")); continue; }
    parts.push(`${d ? d.label : k} ${lo !== null ? "≥ " + lo : ""}${lo !== null && hi !== null ? " and " : ""}${hi !== null ? "≤ " + hi : ""}`);
  }
  return parts.join(" · ");
}
/* does an ion (its descriptor row) already break the body rules? returns the broken rule names */
function brokenBodyRules(desc, c) {
  if (!c || !desc) return [];
  const out = [];
  for (const [k, v] of Object.entries(c)) {
    if (!k.startsWith("body_")) continue;
    if (k === "body_elements") { const bad = (desc.body_elements || []).filter((e) => !v.includes(e)); if (bad.length) out.push(`its body already holds ${bad.join(", ")}`); continue; }
    const d = S.boot.descriptors.find((x) => x.key === k); const [lo, hi] = v; const val = desc[k];
    if (val === undefined) continue;
    if ((hi !== null && val > hi) || (lo !== null && val < lo)) out.push(`it already has ${val} ${d ? d.label.toLowerCase() : k}`);
  }
  return out;
}

function wireRange(id, onChange, unit = "") {
  const r = $("#" + id), o = $("#" + id + "-out");
  const paint = () => { r.style.setProperty("--p", (100 * (r.value - r.min) / (r.max - r.min)) + "%"); o.textContent = r.value + unit; };
  r.addEventListener("input", () => { paint(); if (onChange) onChange(Number(r.value)); });
  paint();
}

/* ---------------- navigation ---------------- */
function show(screen) {
  S.screen = screen;
  $$(".screen").forEach((s) => s.classList.toggle("on", s.id === "screen-" + screen));
  $$("nav button").forEach((b) => b.classList.toggle("on", b.dataset.screen === screen));
  window.scrollTo({ top: 0 });
  try { if (location.hash !== "#" + screen) history.replaceState(null, "", "#" + screen); } catch (_) {}
  const r = RENDER[screen]; if (r) r();
}
document.querySelector("header .brand").addEventListener("click", (e) => { e.preventDefault(); show("start"); });
$$("nav button").forEach((b) => b.addEventListener("click", () => show(b.dataset.screen)));
window.addEventListener("hashchange", () => { const h = location.hash.replace("#", ""); if (RENDER[h] && h !== S.screen) show(h); });

/* ---------------- start ---------------- */
function renderStart() {
  const b = S.boot, c = b.counts, h = b.headline;
  $("#screen-start").innerHTML = `
  <div class="hero"><div class="kicker">What would you like to start from?</div>
  <h1 class="t">Make new ions that keep their family.</h1>
  <p class="lead">An edit can produce a perfectly valid molecule and still destroy the thing that made the ion useful: the charged group that defines its family. IonWeaver finds the charge-bearing core the members of a family share, protects every atom and every bond inside it, and changes only the neutral body around it. Each edit is decided on the parent graph before anything is built, and every product is certified five ways before it is offered. Nothing is repaired, retried or filtered afterwards.</p></div>
  <div class="choices">
    <button class="choice" data-go="families"><span class="num">1</span><strong>A published family</strong><span>Pick one of ${c.families} families of documented ions (${fmt(c.editable_ions)} editable), choose an ion, say how many new ions you want, filter them, look at them, download them.</span></button>
    <button class="choice" data-go="own"><span class="num">2</span><strong>Your own ion</strong><span>Paste one SMILES string. IonWeaver names its family, shows the core it will protect, and offers the other published cores it also matches.</span></button>
    <button class="choice" data-go="sheet"><span class="num">3</span><strong>Your own family</strong><span>Upload a sheet with two columns, <code>SMILES</code> and <code>family</code>. Two ions of a family are enough: the shared core is discovered, then every member can be edited.</span></button>
  </div>
  <div class="card studio-promo"><div class="two"><div><h3>The studio · build an ion by hand</h3><p class="note" style="margin-top:0">Point at an atom or a bond, choose a move, and the engine makes it or names the rule that stops it. Then keep going from the new ion.</p></div><div style="display:grid;place-items:center"><button class="btn" data-go="studio" style="font-size:17px;padding:14px 26px">Open the studio</button></div></div></div>
  <div class="card"><h3>Colour means one thing</h3><div class="roles">${b.roles.map((r) => `<span class="role"><i style="background:${r.color}"></i>${esc(r.label)}</span>`).join("")}</div>
  <p class="note">Families are not given colours of their own: family identity is carried by the blue core itself. Operators are not coloured either, because every operator carries the same guarantee.</p></div>
  <div class="card"><h3>The submitted study, in numbers</h3>
  <div class="tiles">
    <div class="tile"><div class="k">curated ions</div><div class="v">${fmt(c.ions)}</div></div>
    <div class="tile"><div class="k">families</div><div class="v">${c.families}</div></div>
    <div class="tile"><div class="k">protected cores</div><div class="v">${c.cores}</div></div>
    <div class="tile ink"><div class="k">library components</div><div class="v">${c.library_components}</div></div>
    <div class="tile amber"><div class="k">definitions offered</div><div class="v">${fmt(h.definitions_offered)}</div></div>
    <div class="tile green"><div class="k">products certified</div><div class="v">${fmt(h.products_certified_at_the_reported_allowance)}</div></div>
    <div class="tile green"><div class="k">refused by certification</div><div class="v">${h.products_refused_by_certification}</div></div>
    <div class="tile green"><div class="k">validity</div><div class="v">${pct(h.validity_percent, 1)}</div></div>
    <div class="tile ink"><div class="k">novelty vs inventory</div><div class="v">${pct(h.novelty_percent, 1)}</div></div>
  </div><p class="note small">Every value is read from the study's own archive tables (<code>FACTS.json</code>), not typed in. The engine running here is the notebook's, module for module; the Method page shows the hashes.</p></div>
  <div class="two"><div class="card"><h3>Across families</h3><p class="note" style="margin-top:0">Need a set for screening? <a href="#" data-go="batch">Generate a quota of distinct ions per family</a> across every published family, or across the families of your own sheet, with your rules set beforehand, and download one CSV with structures, provenance, the five checks and RDKit descriptors.</p></div>
  <div class="card"><h3>The component library</h3><p class="note" style="margin-top:0">Graft, Swap and Cross-over draw on <a href="#" data-go="library">${c.library_components} components</a> cut from the curated ions themselves, each carried whole with the bond it was attached by. Nothing is invented.</p></div></div>`;
  $$("[data-go]", $("#screen-start")).forEach((el) => el.addEventListener("click", (e) => { e.preventDefault(); show(el.dataset.go); }));
}

/* ---------------- families ---------------- */
function renderFamilies() {
  const root = $("#screen-families");
  if (S.family) { renderFamily(S.family); return; }
  root.innerHTML = `<div class="kicker">Published families</div><h2 class="st">Choose a family</h2>
  <p class="lead">Each family is a group of curated ions sharing one charged core, found by the study's layered common-substructure search. The card shows the family's smallest editable ion with its protected core in blue.</p>
  <div class="filter-row"><span class="chip tog on" data-f="all">All</span><span class="chip tog" data-f="cation">Cations</span><span class="chip tog" data-f="anion">Anions</span><input type="text" id="famSearch" placeholder="search a family…"></div>
  <div class="grid" id="famGrid"></div>`;
  let f = "all", q = "";
  const paint = () => {
    const rows = S.boot.families.filter((x) => (f === "all" || x.ion_type === f) && x.family_label.toLowerCase().includes(q));
    $("#famGrid").innerHTML = rows.map((x, i) => `<div class="gcard" data-id="${esc(x.family_id)}" style="animation-delay:${Math.min(i * 25, 400)}ms">
      <div class="pic">${x.representative_svg || ""}</div>
      <div class="ttl"><span>${esc(x.family_label)}</span>${ionType(x.ion_type)}</div>
      <div class="meta"><span>${x.members} ions</span><span>${x.editable} editable</span><span>${x.cores.length} core${x.cores.length === 1 ? "" : "s"}</span></div></div>`).join("") || `<div class="empty">No family matches.</div>`;
    $$("#famGrid .gcard").forEach((c) => c.addEventListener("click", () => { S.family = c.dataset.id; renderFamily(c.dataset.id); }));
  };
  $$(".filter-row .tog", root).forEach((c) => c.addEventListener("click", () => { $$(".filter-row .tog", root).forEach((x) => x.classList.remove("on")); c.classList.add("on"); f = c.dataset.f; paint(); }));
  $("#famSearch").addEventListener("input", (e) => { q = e.target.value.trim().toLowerCase(); paint(); });
  paint();
}
async function renderFamily(id) {
  const root = $("#screen-families");
  const fam = S.boot.families.find((x) => x.family_id === id);
  root.innerHTML = `<div class="crumbs"><button id="backFams">Families</button><span>${esc(fam.family_label)}</span></div>
  <h2 class="st">${esc(fam.family_label)} ${ionType(fam.ion_type)}</h2>
  <div class="card"><h3>Protected core${fam.cores.length > 1 ? "s" : ""}</h3>
  <table><thead><tr><th>core</th><th>pattern (SMARTS)</th><th class="num">atoms</th><th class="num">bonds</th><th class="num">ions it reaches</th><th class="num">share of family</th></tr></thead>
  <tbody>${fam.cores.map((c) => `<tr><td><b>${c.grammar_id}</b></td><td class="mono">${esc(c.pattern)}</td><td class="num">${c.atoms}</td><td class="num">${c.bonds}</td><td class="num">${c.support}</td><td class="num">${pct(100 * c.support_fraction, 0)}</td></tr>`).join("")}</tbody></table>
  <p class="note small">A family with more than one layer holds ions that share no single core: the best-supported core is G0, the ions it reaches are set aside, and the search repeats on the rest.</p></div>
  ${fam.editable ? `<div class="card"><div class="results-head"><h2 class="st" style="margin:0">New ions from this family</h2><span class="small">no need to pick a parent: ${fam.editable} editable ions, spread from smallest to largest, serve as parents</span></div>` : `<div class="card"><p class="note watch" style="margin:0"><b>No editable ion.</b> No core reached the members of this family in the study, so nothing can be edited here.</p></div><div hidden>`}
    <div class="two" style="margin-top:8px"><div>${rangeKnob("famWant", "How many new ions", 1, 500, 50, 1, "distinct, certified, none equal to a parent or to each other")}${rangeKnob("famParents", "Parents used", 1, Math.max(1, Math.min(40, fam.editable)), Math.min(8, fam.editable), 1, "more parents, more variety")}</div>
    <div>${constraintsPanel("consFam")}<div class="actions"><button class="btn go" id="goFam">Generate</button></div><div id="famProg" hidden><div class="progress"><i></i></div><div class="small pmsg"></div></div></div></div>
    <div id="famOut"></div></div>
  <div class="results-head"><h2 class="st" style="margin:0">Or pick one parent <span class="small">· ${fam.editable} editable of ${fam.members} · click an ion to open it in the workbench</span></h2><div class="actions" style="margin:0"><input type="text" id="ionSearch" placeholder="filter by SMILES…" style="width:220px"></div></div>
  <div class="grid" id="ionGrid"><div class="empty">Drawing ${fam.editable} ions…</div></div>`;
  $("#backFams").addEventListener("click", () => { S.family = null; renderFamilies(); });
  wireRange("famWant", null); wireRange("famParents", null); wireConstraints("consFam");
  if (S.famResult && S.famResult.familyId === id) renderResults($("#famOut"), S.famResult, S.famResult.opts);
  $("#goFam").addEventListener("click", async () => {
    const prog = $("#famProg"), btn = $("#goFam"); prog.hidden = false; btn.disabled = true;
    const per = Number($("#famWant").value);
    const body = { family_ids: [id], per_family: per, ions_per_family: Number($("#famParents").value), constraints: readConstraints("consFam"), svg: per <= 400 };
    try {
      const b = await job("/api/batch", body, (st) => progressBar(prog, st));
      const row = b.families[0] || {};
      S.famResult = { familyId: id, products: b.products, summary: { delivered: b.summary.delivered, requested: per, parents: row.ions_used, new_to_inventory: b.summary.new_to_inventory, elapsed_s: b.summary.elapsed_s }, request: { constraints: body.constraints } };
      S.famResult.opts = { exportName: `ionweaver_${fam.family_label.replace(/\W+/g, "_")}`, tiles: [["green", "distinct certified ions", b.summary.delivered], ["", "parents used", row.ions_used], ["ink", "new to the inventory", b.summary.new_to_inventory], ["", "time", b.summary.elapsed_s + " s"]], note: row.delivered < per ? `<p class="note watch"><b>${row.delivered} of ${per}.</b> These parents offer no further distinct certified ions under these rules${row.parents_breaking_rules ? ` (${row.parents_breaking_rules} of them already break a rule)` : ""}; raise the number of parents.</p>` : "" };
      renderResults($("#famOut"), S.famResult, S.famResult.opts);
    } catch (e) { $("#famOut").innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
    prog.hidden = true; btn.disabled = false;
  });
  try {
    const d = await api(`/api/family/${encodeURIComponent(id)}`);
    if (S.family !== id || S.screen !== "families") return;
    const paint = (q = "") => {
      const ions = d.ions.filter((x) => x.smiles.toLowerCase().includes(q));
      $("#ionGrid").innerHTML = ions.map((x, i) => `<div class="gcard" data-id="${x.inventory_id}" style="animation-delay:${Math.min(i * 12, 500)}ms">
        <div class="pic">${x.svg}</div><div class="smi" title="${esc(x.smiles)}">${esc(x.smiles)}</div>
        <div class="meta"><span>#${x.inventory_id}</span><span>core ${x.grammar_id}</span><span>${x.heavy_atoms} atoms</span><span>${x.body_heavy_atoms} outside the core</span></div></div>`).join("") || `<div class="empty">No ion matches.</div>`;
      $$("#ionGrid .gcard").forEach((c) => c.addEventListener("click", () => openBench({ kind: "published", inventory_id: Number(c.dataset.id) })));
    };
    $("#ionSearch").addEventListener("input", (e) => paint(e.target.value.trim().toLowerCase()));
    paint();
  } catch (e) { $("#ionGrid").innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
}

/* ---------------- own ion ---------------- */
function renderOwn() {
  const root = $("#screen-own");
  if (root.dataset.ready) return;
  root.dataset.ready = "1";
  root.innerHTML = `<div class="kicker">Your own ion</div><h2 class="st">Paste one SMILES string</h2>
  <p class="lead">One cation or one anion, written on its own: one connected structure, no counter-ion, net charge +1 or −1 on a single atom. IonWeaver matches it against the ${S.boot.counts.cores} published cores and protects the most specific one it carries; you may choose another it also matches.</p>
  <div class="card"><div class="field"><label for="ownSmiles">SMILES</label><input type="text" id="ownSmiles" placeholder="e.g. CCCCn1cc[n+](C)c1  or  CCCCCCCC[P+](CCCC)(CCCC)CCCC  or  CC(=O)[O-]" class="mono"></div>
  <div class="actions"><button class="btn" id="ownGo">Identify the family and core</button><span class="small">Try: <a href="#" data-ex="CCCCn1cc[n+](C)c1">1-butyl-3-methylimidazolium</a> · <a href="#" data-ex="CCCCCCCCCCCC[NH3+]">dodecylammonium</a> · <a href="#" data-ex="CCCCC[n+]1ccc(C)cc1">pentyl-picolinium</a> · <a href="#" data-ex="CCCCCCS(=O)(=O)[O-]">hexanesulfonate</a> · <a href="#" data-ex="O=C([O-])c1ccccc1">benzoate</a></span></div>
  <div id="ownOut"></div></div>`;
  $$("[data-ex]", root).forEach((a) => a.addEventListener("click", (e) => { e.preventDefault(); $("#ownSmiles").value = a.dataset.ex; identifyOwn(); }));
  $("#ownGo").addEventListener("click", identifyOwn);
  $("#ownSmiles").addEventListener("keydown", (e) => { if (e.key === "Enter") identifyOwn(); });
}
async function identifyOwn() {
  const out = $("#ownOut"), smiles = $("#ownSmiles").value.trim();
  out.innerHTML = `<p class="note">Reading the structure and matching ${S.boot.counts.cores} cores…</p>`;
  try {
    const v = await api("/api/identify", { smiles });
    if (!v.ok) { out.innerHTML = `<div class="note bad"><b>Not editable as given.</b> ${esc(v.reason)}</div>`; return; }
    const b = v.best;
    out.innerHTML = `<div class="stage" style="margin-top:16px"><div><div class="seedpic">${v.svg}</div>
      <dl class="facts"><dt>canonical SMILES</dt><dd class="mono">${esc(v.canonical_smiles)}</dd><dt>ion type</dt><dd>${v.ion_type}, charge ${v.charge > 0 ? "+1" : "−1"}</dd>
      <dt>best match</dt><dd><b>${esc(b.family_label)}</b> · core ${b.grammar_id} (${b.atoms} atoms, carried by ${b.support} published ions)</dd><dt>core pattern</dt><dd class="mono">${esc(b.pattern)}</dd>
      ${v.published_id !== null ? `<dt>in the inventory</dt><dd>yes, reference ion ${v.published_id}; edited under the core the study assigned it</dd>` : ""}</dl></div>
      <div><div class="card" style="margin:0"><h3>Protect this core</h3><p class="note" style="margin-top:0">Atoms Q0, Q1, … in blue are held through every edit. Everything else may change, including a group you may think of as part of the family (an S=O on a sulfoxonium matched as a sulfonium, for instance). If the defining group is larger than this core, <a href="#" data-go="sheet">upload a sheet of the family</a> and its own core is discovered.</p>
      <div class="actions"><button class="btn" id="useBest">Edit under ${esc(b.family_label)} ${b.grammar_id}</button></div>
      ${v.alternatives.length ? `<h3 style="margin-top:18px">Other cores it also carries</h3><div class="grid" style="grid-template-columns:repeat(auto-fill,minmax(170px,1fr))">${v.alternatives.map((a, i) => `<div class="gcard" data-alt="${i}"><div class="pic" style="min-height:110px">${a.svg}</div><div class="ttl" style="font-size:15px"><span>${esc(a.family_label)} ${a.grammar_id}</span></div><div class="meta"><span>${a.atoms} atoms</span><span>${a.support} ions</span></div></div>`).join("")}</div>` : `<p class="note small">No other published core matches this ion.</p>`}</div></div></div>`;
    const go = (m) => openBench({ kind: "own", smiles: v.canonical_smiles, family_id: m.family_id, family_label: m.family_label, grammar_id: m.grammar_id, pattern: m.pattern, core_atoms: m.core_atoms });
    $$("[data-go]", out).forEach((a) => a.addEventListener("click", (e) => { e.preventDefault(); show(a.dataset.go); }));
    $("#useBest").addEventListener("click", () => go(b));
    $$("[data-alt]", out).forEach((c) => c.addEventListener("click", () => go(v.alternatives[Number(c.dataset.alt)])));
  } catch (e) { out.innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
}

/* ---------------- own family (sheet) ---------------- */
function renderSheet() {
  const root = $("#screen-sheet");
  if (root.dataset.ready) { return; }
  root.dataset.ready = "1";
  root.innerHTML = `<div class="kicker">Your own family</div><h2 class="st">Upload a sheet of SMILES and family labels</h2>
  <p class="lead">Two columns, <code>SMILES</code> and <code>family</code>; an <code>ion type</code> or <code>charge</code> column is used when present. Two surviving ions make a family. Each family's shared core is discovered with the study's own search, and families nobody has published get a core of their own.</p>
  <div class="two"><div class="card"><h3>File</h3><div class="drop" id="drop">Drop a .csv, .tsv, .txt or .xlsx here, or click to choose<br><span class="small">up to 3,000 rows · nothing is stored beyond this session</span><input type="file" id="fileIn" accept=".csv,.tsv,.txt,.xlsx" hidden></div>
    <div class="actions"><a class="btn sec sm" href="${API}/api/template.csv">Download the template</a><button class="btn sec sm" id="sheetExample">Try the template as a sheet</button><button class="btn sec sm" id="sheetNovel">Try families nobody published</button></div>
    <div class="small" style="margin-top:8px">The second example holds pyrylium, sulfoxonium, dithiocarbamate and iodonium families, none of them in the published inventory, plus two cases that are reported rather than edited: a family of one ion, and plain thiolates whose shared part is too small for a core.</div></div>
  <div class="card"><h3>Or paste rows</h3><textarea id="sheetText" placeholder="SMILES,family&#10;CCCCn1cc[n+](C)c1,imidazolium&#10;CCCCCCn1cc[n+](C)c1,imidazolium"></textarea><div class="actions"><button class="btn" id="sheetPaste">Curate and discover cores</button></div>
    <details class="guide"><summary>What makes a sheet work</summary><ul class="small">
      <li><b>One ion per row, written on its own.</b> No counter-ion in the SMILES; a salt is reported as disconnected.</li>
      <li><b>Give a family at least two ions,</b> and let them differ in their side chains. The core is the largest structure every member shares; two near-identical members make it almost the whole ion, which leaves little to edit.</li>
      <li><b>Members must share more than the charged atom and one neighbour.</b> A core has at least three atoms and two bonds, so plain thiolates (C–S⁻) are reported as uncovered.</li>
      <li><b>Labels are yours.</b> A label is only a grouping; a published name does not import the published core. Cations and anions under one label are two families.</li>
      <li><b>Net charge must be +1 or −1 on one atom.</b> Divalent ions, radicals and metals are set aside, each with its reason, in the row-by-row report.</li>
    </ul></details></div></div>
  <div id="sheetProg" hidden><div class="progress"><i></i></div><div class="small pmsg"></div></div>
  <div id="sheetOut"></div>`;
  const drop = $("#drop"), fin = $("#fileIn");
  drop.addEventListener("click", () => fin.click());
  drop.addEventListener("dragover", (e) => { e.preventDefault(); drop.classList.add("over"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("over"));
  drop.addEventListener("drop", (e) => { e.preventDefault(); drop.classList.remove("over"); if (e.dataTransfer.files[0]) sendFile(e.dataTransfer.files[0]); });
  fin.addEventListener("change", () => { if (fin.files[0]) sendFile(fin.files[0]); });
  $("#sheetExample").addEventListener("click", () => runSheet({ example: true }));
  $("#sheetNovel").addEventListener("click", () => runSheet({ example: "novel" }));
  $("#sheetPaste").addEventListener("click", () => { const t = $("#sheetText").value.trim(); if (!t) { toast("Paste some rows first"); return; } runSheet({ text: t, filename: "pasted.csv" }); });
  function sendFile(f) { const rd = new FileReader(); rd.onload = () => runSheet({ filename: f.name, content_base64: rd.result.split(",")[1] }); rd.readAsDataURL(f); }
}
async function runSheet(body) {
  const prog = $("#sheetProg"), out = $("#sheetOut");
  prog.hidden = false; out.innerHTML = "";
  try {
    const d = await job("/api/sheet", body, (st) => progressBar(prog, st));
    S.sheet = d; S.batchSheet = null; S.batch = null; prog.hidden = true; paintSheet(d);
  } catch (e) { prog.hidden = true; out.innerHTML = `<div class="note bad"><b>The sheet could not be used.</b> ${esc(e.message)}</div>`; }
}
function paintSheet(d) {
  const out = $("#sheetOut");
  const statuses = Object.entries(d.status_counts).map(([k, v]) => `<span class="chip ${k === "retained" ? "green" : "red"}">${esc(v.text)} · ${v.n}</span>`).join("");
  out.innerHTML = `<div class="card"><h3>${esc(d.filename)} · curated in ${d.elapsed_s} s</h3>
  <div class="tiles"><div class="tile"><div class="k">rows read</div><div class="v">${d.rows}</div></div><div class="tile green"><div class="k">ions kept</div><div class="v">${d.retained}</div></div><div class="tile"><div class="k">families</div><div class="v">${d.families.length}</div></div><div class="tile"><div class="k">editable ions (a core was found)</div><div class="v">${d.editable}</div></div><div class="tile ink"><div class="k">components in the sheet</div><div class="v">${d.library_components}</div></div><div class="tile ink"><div class="k">of them new to the published library</div><div class="v">${d.library_new_to_published}</div></div></div>
  <div class="status-list">${statuses}</div>
  <details style="margin-top:10px"><summary class="small" style="cursor:pointer">Row by row</summary><div class="audit"><table><thead><tr><th>row</th><th>SMILES</th><th>family</th><th>type</th><th>outcome</th></tr></thead><tbody>${d.audit.map((r) => `<tr><td>${r.row}</td><td class="mono">${esc(r.smiles)}</td><td>${esc(r.family)}</td><td>${esc(r.ion_type)}</td><td>${r.status === "retained" ? '<span class="chip green" style="padding:2px 9px;font-size:12px">kept</span>' : `<span class="chip red" style="padding:2px 9px;font-size:12px">${esc(r.text)}</span>`}</td></tr>`).join("")}</tbody></table></div></details>
  <div class="actions">${d.editable ? `<button class="btn go" id="sheetBatch">Generate across this sheet's families</button>` : `<span class="small">No editable ion in this sheet, so nothing to generate.</span>`}<span class="small">Click an ion below to open it in the workbench.</span></div></div>
  ${d.families.map((f) => `<div class="card"><div class="results-head"><h2 class="st" style="margin:0">${esc(f.family_label)} ${ionType(f.ion_type)} <span class="small">· ${f.members} kept, ${f.editable} editable</span></h2></div>
    ${f.cores.length ? f.cores.map((c) => `<div class="stage" style="margin-top:12px"><div class="seedpic" style="min-height:160px">${c.svg || ""}</div><div>
      <dl class="facts" style="margin-top:0"><dt>core</dt><dd><b>${c.grammar_id}</b> · ${c.atoms} atoms, ${c.bonds} bonds</dd><dt>pattern</dt><dd class="mono">${esc(c.pattern)}</dd><dt>reaches</dt><dd>${c.support} of ${f.members} ions (${pct(100 * c.support_fraction, 0)})</dd>
      <dt>published core?</dt><dd>${c.published_equivalent ? `yes · the same pattern as <b>${esc(c.published_equivalent.family_label)} ${c.published_equivalent.grammar_id}</b>` : "no · a core of your own"}</dd></dl></div></div>`).join("")
      : `<p class="note watch"><b>No core was found.</b> ${f.members < 2 ? "A family needs at least two kept ions." : "No common substructure met the floors (3 atoms, 2 bonds, the charged atom inside, every charge on it)."}</p>`}
    ${f.uncovered.length ? `<p class="note small">Not reached by any core: ${f.uncovered.map((u) => `<code>${esc(u.smiles)}</code> (${esc(u.reason)})`).join(", ")}</p>` : ""}
    <div class="grid" style="margin-top:12px">${d.ions.filter((x) => x.family_id === f.family_id).map((x) => `<div class="gcard" data-id="${x.inventory_id}"><div class="pic">${x.svg}</div><div class="smi">${esc(x.canonical_smiles)}</div><div class="meta"><span>core ${x.grammar_id}</span><span>${x.heavy_atoms} atoms</span><span>${x.body_heavy_atoms} outside the core</span></div></div>`).join("")}</div></div>`).join("")}`;
  $$(".gcard[data-id]", out).forEach((c) => c.addEventListener("click", () => openBench({ kind: "sheet", sheet_id: d.sheet_id, inventory_id: Number(c.dataset.id) })));
  if ($("#sheetBatch")) $("#sheetBatch").addEventListener("click", () => { S.batchSheet = d; show("batch"); });
}

/* ---------------- workbench ---------------- */
async function openBench(start) {
  S.start = start; S.seed = null; S.gen = null; S.chain = null; S.filters = null;
  $("#navBench").disabled = false;
  show("bench");
}
function startName() {
  const s = S.start; if (!s) return "";
  if (s.kind === "published") return `reference ion ${s.inventory_id}`;
  if (s.kind === "own") return `your ion under ${s.family_label} ${s.grammar_id}`;
  if (s.kind === "state") return `the ion you built (${s.family_label} ${s.grammar_id})`;
  return `sheet ion ${s.inventory_id}`;
}
async function renderBench() {
  const root = $("#screen-bench");
  if (!S.start) { root.innerHTML = `<div class="kicker">Workbench</div><h2 class="st">Choose a parent first</h2><div class="card"><div class="actions" style="margin-top:0"><button class="btn sec" data-go="families">A published family</button><button class="btn sec" data-go="own">Your own ion</button><button class="btn sec" data-go="sheet">Your own family</button></div></div>`; $$("[data-go]", root).forEach((b) => b.addEventListener("click", () => show(b.dataset.go))); return; }
  if (!S.seed) {
    const mine = S.start;
    root.innerHTML = `<div class="kicker">Workbench</div><h2 class="st">Preparing ${esc(startName())}…</h2><div class="progress busy"><i style="width:30%"></i></div><p class="small">deciding every definition the eight operators offer</p>`;
    try { const seed = await api("/api/start", mine); if (S.start !== mine) return; S.seed = seed; } catch (e) { if (S.start === mine) root.innerHTML = `<div class="note bad">${esc(e.message)}</div>`; return; }
  }
  const s = S.seed, ops = S.boot.operators;
  root.innerHTML = `<div class="kicker">Workbench · ${esc(s.workspace === "published" ? "published inventory" : "published inventory + your sheet")}</div>
  <h2 class="st">${esc(s.family_label || s.family_id)} ${ionType(s.ion_type)} <span class="small">· core ${s.grammar_id}${s.inventory_id >= 0 ? " · #" + s.inventory_id : ""}</span></h2>
  <div class="stage">
    <div class="card"><h3>Parent</h3><div class="small" style="margin:-6px 0 8px">The core in blue is held through every edit; everything else may change.</div><div class="seedpic">${s.svg}</div>
      <dl class="facts"><dt>SMILES</dt><dd class="mono">${esc(s.smiles)}</dd><dt>core pattern</dt><dd class="mono">${esc(s.core_pattern)}</dd><dt>size</dt><dd>${s.heavy_atoms} atoms (H not counted), ${s.body_heavy_atoms} outside the core${s.ring_shared_atoms.length ? `, ${s.ring_shared_atoms.length} held in a ring the core covers in part` : ""}</dd>
      <dt>molecular weight</dt><dd>${fmt(s.descriptors.mw, 1)} g/mol · ${esc(s.descriptors.formula)}</dd></dl>
      <h3 style="margin-top:16px">Parts</h3><div class="small" style="margin:-6px 0 8px">The core, and the side groups an edit may swap, cut or splice.</div>
      <div class="parts"><div class="part-card core"><div class="pic">${s.core_svg || ""}</div><div class="k">protected core · ${s.core_atoms.length} atoms</div><div class="small mono">${esc(s.core_pattern)}</div></div>
      ${s.side_groups.length ? s.side_groups.map((g) => `<div class="part-card"><div class="pic">${g.svg || ""}</div><div class="k">side group at Q${s.core_atoms.indexOf(g.anchor_core_atom)} · ${g.heavy_atoms} atom${g.heavy_atoms === 1 ? "" : "s"}</div><div class="small mono">${esc(g.key)}</div></div>`).join("") : `<div class="part-card"><div class="empty" style="padding:16px">no side group: every atom is core${s.ring_shared_atoms.length ? " or held in a ring the core covers in part" : ""}; Insert and Graft can still add at a free valence</div></div>`}</div>
      <h3 style="margin-top:16px">Operators · ${fmt(s.total_capacity)} accepted definitions</h3><div class="small" style="margin:-6px 0 8px">What each operator offers on this parent; click one to leave it out.</div>
      <div class="ops" id="ops">${ops.map((o) => { const n = s.capacity[o.key] || 0; return `<div class="op ${n ? "on" : "none"}" data-op="${o.key}" title="${esc(o.text)}"><div class="k"><span>${o.key}</span><b>${fmt(n)}</b></div><div class="d">${(s.modes[o.key] || []).map((m) => `${esc(m.label)} ${m.count}`).join(" · ") || "no site on this ion"}</div></div>`; }).join("")}</div></div>
    <div><div class="card"><div class="chips" style="margin:0 0 12px"><span class="chip tog ${S.bench === "generate" ? "on" : ""}" data-mode="generate">New ions, one edit each</span><span class="chip tog ${S.bench === "chain" ? "on" : ""}" data-mode="chain">Edit again and again</span></div>
      <div id="modeGenerate" ${S.bench === "generate" ? "" : "hidden"}>
        ${rangeKnob("want", "How many new ions", 1, 500, 60, 1, "distinct, certified, none identical to the parent or to each other")}
        <div class="field"><label>Balance</label><select id="balance"><option value="mix">Interleave the operators (a mix)</option><option value="per">Cap per operator</option></select></div>
        <div id="perWrap" hidden>${rangeKnob("perOp", "Cap per operator", 1, 100, 10, 1, "")}</div>
        ${constraintsPanel("consGen")}
        <div class="actions"><button class="btn go" id="goGen">Generate</button><button class="btn sec" id="benchStudio">Build by hand instead</button></div>
        <p class="small" style="margin:10px 0 0">Candidates come in the study's own examination order; "more" walks further down it.</p></div>
      <div id="modeChain" ${S.bench === "chain" ? "" : "hidden"}>
        <div class="field"><label>Selection rule</label><select id="rule">${S.boot.rules.map((r) => `<option value="${r.key}">${esc(r.label)}</option>`).join("")}</select><div class="hint small" id="ruleHint">${esc(S.boot.rules[0].text)}</div></div>
        ${rangeKnob("depth", "Rounds", 1, 12, 6, 1, "each selected product becomes the next starting structure; the original core stays the reference")}
        ${rangeKnob("allow", "Candidates examined per round", 8, 200, 50, 1, "in the study's own order, across the chosen operators")}
        <div class="actions"><button class="btn go" id="goChain">Run the chain</button></div></div>
      <div id="benchProg" hidden><div class="progress"><i></i></div><div class="small pmsg"></div></div></div></div>
  </div>
  <div id="benchOut"></div>`;
  wireRange("want", null); wireRange("perOp", null); wireRange("depth", null); wireRange("allow", null); wireConstraints("consGen");
  $("#balance").addEventListener("change", (e) => { $("#perWrap").hidden = e.target.value !== "per"; });
  $("#rule").addEventListener("change", (e) => { $("#ruleHint").textContent = (S.boot.rules.find((r) => r.key === e.target.value) || {}).text || ""; });
  $$("#ops .op").forEach((o) => o.addEventListener("click", () => { if (!o.classList.contains("none")) o.classList.toggle("on"); }));
  $$("[data-mode]", root).forEach((c) => c.addEventListener("click", () => { S.bench = c.dataset.mode; $$("[data-mode]", root).forEach((x) => x.classList.toggle("on", x === c)); $("#modeGenerate").hidden = S.bench !== "generate"; $("#modeChain").hidden = S.bench !== "chain"; }));
  $("#goGen").addEventListener("click", runGenerate);
  $("#goChain").addEventListener("click", runChain);
  $("#benchStudio").addEventListener("click", () => openStudio(S.start));
  if (S.gen) paintGenerate(); else if (S.chain) paintChain();
}
const chosenOps = () => $$("#ops .op.on").map((o) => o.dataset.op);
async function runGenerate() {
  const ops = chosenOps(); if (!ops.length) { toast("Include at least one operator"); return; }
  const prog = $("#benchProg"), btn = $("#goGen"); prog.hidden = false; btn.disabled = true;
  const mine = S.start;
  const body = { start: mine, want: Number($("#want").value), operators: ops, per_operator: $("#balance").value === "per" ? Number($("#perOp").value) : null, constraints: readConstraints("consGen") };
  const broken = brokenBodyRules(S.seed && S.seed.descriptors, body.constraints);
  try { const g = await job("/api/generate", body, (st) => progressBar(prog, st)); if (S.start !== mine) return; S.gen = g; S.gen.request = body; S.gen.broken = broken; S.chain = null; paintGenerate(); }
  catch (e) { if (S.start === mine) $("#benchOut").innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
  if (S.start === mine) { prog.hidden = true; btn.disabled = false; }
}
async function runMore() {
  const g = S.gen, prog = $("#benchProg"), btn = $("#moreBtn"); if (!g) return;
  prog.hidden = false; btn.disabled = true;
  const body = { ...g.request, want: Number($("#want").value), cursor: g.summary.cursor, exclude: g.products.map((p) => p.product_smiles) };
  try {
    const more = await job("/api/generate", body, (st) => progressBar(prog, st));
    if (S.gen !== g) return;
    const s = g.summary, m = more.summary;
    g.products = g.products.concat(more.products);
    g.summary = { ...m, delivered: g.products.length, requested: s.requested + m.requested, examinations: s.examinations + m.examinations,
      refused_by_certification: s.refused_by_certification + m.refused_by_certification, repeats: s.repeats + m.repeats,
      identical_to_parent: s.identical_to_parent + m.identical_to_parent, outside_constraints: (s.outside_constraints || 0) + (m.outside_constraints || 0),
      new_to_inventory: s.new_to_inventory + m.new_to_inventory, elapsed_s: Math.round(100 * (s.elapsed_s + m.elapsed_s)) / 100, cursor: m.cursor, exhausted: m.exhausted };
    g.filters = null; paintGenerate();
    if (!more.products.length) toast("No further distinct certified ion is available from one edit");
  } catch (e) { toast(e.message); }
  prog.hidden = true; btn.disabled = false;
}
async function runChain() {
  const ops = chosenOps(); if (!ops.length) { toast("Include at least one operator"); return; }
  const prog = $("#benchProg"), btn = $("#goChain"); prog.hidden = false; btn.disabled = true;
  const mine = S.start;
  const body = { start: mine, rule: $("#rule").value, depth: Number($("#depth").value), allowance: Number($("#allow").value), operators: ops };
  try { const c = await job("/api/chain", body, (st) => progressBar(prog, st)); if (S.start !== mine) return; S.chain = c; S.gen = null; paintChain(); }
  catch (e) { if (S.start === mine) $("#benchOut").innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
  if (S.start === mine) { prog.hidden = true; btn.disabled = false; }
}

/* ----- results: filters, sort, cards ----- */
function buildFilters(products) {
  const f = { ops: new Set(products.map((p) => p.operator)), newOnly: false, ranges: {} };
  for (const d of S.boot.descriptors) {
    const vals = products.map((p) => p.descriptors[d.key]).filter((v) => v !== undefined && v !== null);
    if (!vals.length) continue;
    const lo = Math.min(...vals), hi = Math.max(...vals);
    if (lo === hi) continue;
    f.ranges[d.key] = { lo, hi, min: lo, max: hi, step: d.decimals ? Math.pow(10, -d.decimals) : 1 };
  }
  const sims = products.map((p) => p.similarity_to_parent);
  if (sims.length && Math.min(...sims) !== Math.max(...sims)) f.ranges.similarity_to_parent = { lo: Math.min(...sims), hi: Math.max(...sims), min: Math.min(...sims), max: Math.max(...sims), step: 0.01 };
  return f;
}
function applyFilters(products, f) {
  return products.filter((p) => f.ops.has(p.operator) && (!f.newOnly || p.new_to_inventory) && Object.entries(f.ranges).every(([k, r]) => {
    const v = k === "similarity_to_parent" ? p.similarity_to_parent : p.descriptors[k];
    return v >= r.min - 1e-9 && v <= r.max + 1e-9;
  }));
}
function sortProducts(products) {
  const [key, dir] = S.sort.split("_").length > 2 ? [S.sort.slice(0, S.sort.lastIndexOf("_")), S.sort.slice(S.sort.lastIndexOf("_") + 1)] : S.sort.split("_");
  const val = (p) => key === "similarity" ? p.similarity_to_parent : key === "operator" ? p.operator : key === "order" ? 0 : p.descriptors[key];
  const out = [...products];
  if (key === "order") return out;
  out.sort((a, b) => (val(a) > val(b) ? 1 : val(a) < val(b) ? -1 : 0) * (dir === "desc" ? -1 : 1));
  return out;
}
function productCard(p, i) {
  return `<div class="pcard" data-i="${i}" style="animation-delay:${Math.min(i * 18, 600)}ms">
    <div class="top">${opchip(p.operator)}<span class="sim" title="Tanimoto similarity to the parent, Morgan radius 2">${p.similarity_to_parent.toFixed(2)}</span></div>
    <div class="pic">${p.svg || ""}</div>
    <div class="what">${esc(p.description)}</div>
    <div class="desc"><span>${esc(p.descriptors.formula)}</span><span>${fmt(p.descriptors.mw, 1)} g/mol</span><span>HBA ${p.descriptors.hba} · HBD ${p.descriptors.hbd}</span>${p.new_to_inventory ? '<span class="newchip" title="not in the reference inventory">new</span>' : '<span class="small" title="already in the reference inventory">known</span>'}</div>
    <div class="smi" title="${esc(p.product_smiles)}">${esc(p.product_smiles)}</div></div>`;
}
function paintGenerate() {
  const g = S.gen, sm = g.summary;
  const tiles = [["green", "distinct certified ions", sm.delivered], ["", "candidates examined", sm.examinations], ["ink", "refused by certification", sm.refused_by_certification], ["amber", "repeats set aside", sm.repeats + sm.identical_to_parent]];
  if (g.request && g.request.constraints) tiles.push(["amber", "outside your rules", sm.outside_constraints || 0]);
  tiles.push(["ink", "new to the inventory", sm.new_to_inventory], ["", "time", sm.elapsed_s + " s"]);
  let note = "";
  if (g.broken && g.broken.length) note += `<p class="note watch"><b>The parent itself breaks a rule:</b> ${esc(g.broken.join("; "))}. Only edits that remove those atoms can comply.</p>`;
  if (sm.exhausted) note += `<p class="note watch"><b>Every definition was examined.</b> This parent offers ${sm.delivered} distinct certified ions under the chosen operators${g.request && g.request.constraints ? " and rules" : ""}; a larger number is not available from one edit.</p>`;
  else if (sm.stopped_on_budget) note += `<p class="note watch"><b>Stopped after ${fmt(sm.examinations)} examinations</b> with ${sm.delivered} ions that obey the rules. "More" keeps looking from where it stopped.</p>`;
  if (!sm.exhausted) note += `<div class="actions"><button class="btn sec" id="moreBtn">${$("#want") ? Number($("#want").value) : 60} more</button><button class="btn sec" id="toStudio">Open the parent in the studio</button></div>`;
  renderResults($("#benchOut"), g, { exportName: `ionweaver_${(g.start.family_label || "ions").replace(/\W+/g, "_")}`, tiles, note, onReset: paintGenerate });
  if ($("#moreBtn")) $("#moreBtn").addEventListener("click", runMore);
  if ($("#toStudio")) $("#toStudio").addEventListener("click", () => openStudio(S.start));
}
/* one results block: summary tiles, filters, sort, cards, downloads. g = {products, summary, request} */
function renderResults(out, g, opts) {
  const all = g.products;
  if (!g.filters) g.filters = buildFilters(all);
  const f = g.filters;
  out.innerHTML = `<div class="card"><h3>What came back</h3>
  <div class="tiles">${opts.tiles.map(([cls, k, v]) => `<div class="tile ${cls}"><div class="k">${k}</div><div class="v ${String(v).endsWith(" s") ? "small" : ""}">${v}</div></div>`).join("")}</div>
  ${g.request && g.request.constraints ? `<p class="note"><b>Rules in force:</b> ${esc(constraintsText(g.request.constraints))}.</p>` : ""}
  ${opts.note || ""}
  <p class="note small">Every ion shown passed all five checks; click a card for its certificate.</p></div>
  <div class="card"><div class="results-head"><div class="small"><b><span class="fcount">${all.length}</span> of ${all.length}</b> shown</div><div class="actions" style="margin:0"><label class="chip tog ${f.newOnly ? "on" : ""} newOnly" title="not in the reference inventory">new only</label><select class="sortSel" style="width:auto"><option value="order">in the examined order</option><option value="similarity_desc">most similar first</option><option value="similarity_asc">least similar first</option><option value="mw_asc">lightest first</option><option value="mw_desc">heaviest first</option><option value="logp_asc">lowest log P first</option><option value="logp_desc">highest log P first</option><option value="operator_asc">by operator</option></select></div></div>
  <div class="chips fops">${S.boot.operators.map((o) => `<span class="chip tog ${f.ops.has(o.key) ? "on" : ""}" data-op="${o.key}">${o.key} <span class="lvl">${all.filter((p) => p.operator === o.key).length}</span></span>`).join("")}</div>
  <div class="filters franges">${Object.entries(f.ranges).map(([k, r]) => { const d = S.boot.descriptors.find((x) => x.key === k) || { label: "Similarity to parent", decimals: 2, unit: "" }; return `<div class="f" data-k="${k}"><label><span>${esc(d.label)}${d.unit ? " (" + d.unit + ")" : ""}</span><b><span class="lo">${fmt(r.min, d.decimals)}</span> – <span class="hi">${fmt(r.max, d.decimals)}</span></b></label>
    <div class="dual"><div class="track"></div><div class="fill"></div><input type="range" class="rlo" min="${r.lo}" max="${r.hi}" step="${r.step}" value="${r.min}"><input type="range" class="rhi" min="${r.lo}" max="${r.hi}" step="${r.step}" value="${r.max}"></div></div>`; }).join("")}</div>
  <div class="actions"><span class="small">Download the shown set</span><button class="btn sm expCsv">CSV</button><button class="btn sec sm expSdf">SDF</button><button class="btn sec sm expSmi">SMILES</button><button class="btn sec sm resetF" style="margin-left:auto">Reset filters</button></div></div>
  <div class="pgrid"></div>`;
  $(".sortSel", out).value = S.sort;
  const paint = () => {
    const shown = sortProducts(applyFilters(all, f));
    $(".fcount", out).textContent = shown.length;
    $(".pgrid", out).innerHTML = shown.map((p) => productCard(p, all.indexOf(p))).join("") || `<div class="empty" style="grid-column:1/-1">Nothing passes these filters.</div>`;
    $$(".pgrid .pcard", out).forEach((c) => c.addEventListener("click", () => openProduct(all[Number(c.dataset.i)])));
    const exp = (fmt_) => exportProducts(shown, fmt_, opts.exportName || "ionweaver");
    $(".expCsv", out).onclick = () => exp("csv"); $(".expSdf", out).onclick = () => exp("sdf"); $(".expSmi", out).onclick = () => exp("smi");
  };
  $(".sortSel", out).addEventListener("change", (e) => { S.sort = e.target.value; paint(); });
  $(".newOnly", out).addEventListener("click", (e) => { f.newOnly = !f.newOnly; e.currentTarget.classList.toggle("on", f.newOnly); paint(); });
  $$(".fops .tog", out).forEach((c) => c.addEventListener("click", () => { const k = c.dataset.op; if (f.ops.has(k)) f.ops.delete(k); else f.ops.add(k); c.classList.toggle("on", f.ops.has(k)); paint(); }));
  $$(".franges .f", out).forEach((box) => {
    const k = box.dataset.k, r = f.ranges[k], lo = $(".rlo", box), hi = $(".rhi", box), fill = $(".fill", box);
    const d = S.boot.descriptors.find((x) => x.key === k) || { decimals: 2 };
    const sync = () => { let a = Number(lo.value), b = Number(hi.value); if (a > b) [a, b] = [b, a]; if (a <= r.lo + r.step) a = r.lo; if (b >= r.hi - r.step) b = r.hi; r.min = a; r.max = b; $(".lo", box).textContent = fmt(a, d.decimals); $(".hi", box).textContent = fmt(b, d.decimals); const span = r.hi - r.lo || 1; fill.style.left = (100 * (a - r.lo) / span) + "%"; fill.style.width = (100 * (b - a) / span) + "%"; };
    lo.addEventListener("input", () => { sync(); paint(); }); hi.addEventListener("input", () => { sync(); paint(); }); sync();
  });
  $(".resetF", out).addEventListener("click", () => { g.filters = null; if (opts.onReset) opts.onReset(); else renderResults(out, g, opts); });
  paint();
}
async function openProduct(p) {
  const dlg = $("#dlg"), body = $("#dlgBody");
  const checks = S.boot.requirements.map((r) => `<div class="check"><i class="${p.checks[r.key] ? "" : "no"}">${p.checks[r.key] ? "✓" : "✕"}</i><div><b>${esc(r.short)}.</b> ${esc(r.text)}</div></div>`).join("");
  const title = p.edit ? p.edit.charAt(0).toUpperCase() + p.edit.slice(1) : p.operator;
  body.innerHTML = `<button class="dlg-close" id="dlgClose">×</button><div class="kicker">${esc(p.operator)}${p.round ? ` · round ${p.round}` : ""} · certified, five of five</div><h2>${esc(title)}</h2><p class="note" style="margin:0 0 12px">${esc(p.description)}. <span class="small">Atom numbers are the small grey labels on the parent below.</span></p>
  <div class="pairs"><div><div class="small">parent · removed in rose, touched in amber</div><div class="pic" id="parentPic"><span class="small">drawing…</span></div></div><div class="arrow">→</div><div><div class="small">product · added or changed in amber</div><div class="pic">${p.svg || ""}</div></div></div>
  <div class="two" style="margin-top:14px"><div><h3 class="kicker">Certification</h3><div class="checks">${checks}</div><p class="small" style="margin-top:8px">An independent judge re-derived the same verdict by searching the product for the core: ${p.judge_agrees ? "it agrees." : "it disagrees (recorded)."}</p></div>
  <div><h3 class="kicker">Provenance</h3><dl class="kv"><dt>parent</dt><dd class="mono">${esc(p.parent_smiles)}</dd><dt>product</dt><dd class="mono">${esc(p.product_smiles)}</dd><dt>mapped product</dt><dd class="mono small">${esc(p.mapped_smiles)}</dd><dt>site</dt><dd>${esc(p.site)}</dd><dt>graph change</dt><dd class="small">${esc(p.signature)}</dd><dt>definition</dt><dd>${esc(p.definition_id)} · examined ${p.order_index ? "at position " + p.order_index : ""}</dd><dt>similarity to parent</dt><dd>${p.similarity_to_parent.toFixed(2)} (Tanimoto, Morgan r=2, 2048 bits)</dd><dt>size change</dt><dd>${p.heavy_atom_change === 0 ? "none" : (p.heavy_atom_change > 0 ? "+" : "") + p.heavy_atom_change + (Math.abs(p.heavy_atom_change) === 1 ? " atom" : " atoms")}</dd><dt>in the inventory</dt><dd>${p.new_to_inventory ? "no · new" : "yes · a known ion"}</dd></dl></div></div>
  <h3 class="kicker" style="margin-top:14px">RDKit descriptors</h3><div class="desc-grid">${S.boot.descriptors.filter((d) => !d.body).map((d) => `<div title="${esc(d.how)}">${esc(d.label)}<b>${fmt(p.descriptors[d.key], d.decimals)}${d.unit ? " " + d.unit : ""}</b></div>`).join("")}</div>
  <h3 class="kicker" style="margin-top:12px">Outside the core</h3><div class="desc-grid"><div>Elements<b>${esc((p.descriptors.body_elements || []).join(", ") || "none")}</b></div>${S.boot.descriptors.filter((d) => d.body).map((d) => `<div title="${esc(d.how)}">${esc(d.label.replace(" outside the core", ""))}<b>${fmt(p.descriptors[d.key], 0)}</b></div>`).join("")}</div>
  <div class="actions"><button class="btn sm" id="copySmi">Copy SMILES</button><button class="btn sec sm" id="useAsStart" title="replaces the current results">Use as the new parent</button><button class="btn sec sm" id="useInStudio">Open in the studio</button></div>`;
  dlg.showModal();
  $("#dlgClose").onclick = () => dlg.close();
  $("#copySmi").onclick = () => { navigator.clipboard.writeText(p.product_smiles); toast("SMILES copied"); };
  const asStart = { kind: "state", smiles: p.product_smiles, core_atoms: p.product_core_atoms || null, family_id: p.family_id, family_label: p.family_label, grammar_id: p.grammar_id, pattern: p.core_pattern, workspace: p.workspace || "published", inventory_id: -1 };
  $("#useAsStart").onclick = () => { dlg.close(); openBench(asStart); };
  $("#useInStudio").onclick = () => { dlg.close(); openStudio(asStart); };
  try {
    const parentCore = p.parent_core_atoms || (S.seed && S.seed.smiles === p.parent_smiles ? S.seed.core_atoms : []);
    const d = await api("/api/depict", { smiles: p.parent_smiles, core_atoms: parentCore, removed_atoms: p.removed_parent_atoms, added_atoms: (p.changed_parent_atoms || []).filter((i) => !parentCore.includes(i)), number_body: true, width: 420, height: 300 });
    $("#parentPic").innerHTML = d.svg;
  } catch (_) { $("#parentPic").innerHTML = `<span class="small">could not draw the parent</span>`; }
}
$("#dlg").addEventListener("click", (e) => { if (e.target === $("#dlg")) $("#dlg").close(); });

function paintChain() {
  const c = S.chain, out = $("#benchOut");
  const moved = c.steps.filter((s) => s.moved);
  out.innerHTML = `<div class="card"><h3>${esc(S.boot.rules.find((r) => r.key === c.rule).label)} · ${c.depth} rounds asked, ${c.summary.rounds_moved} taken</h3>
  <div class="tiles"><div class="tile green"><div class="k">rounds moved</div><div class="v">${c.summary.rounds_moved}</div></div><div class="tile"><div class="k">final size, atoms</div><div class="v">${c.summary.final_heavy_atoms}</div></div><div class="tile amber"><div class="k">similarity, first to last</div><div class="v">${c.summary.similarity_first_to_last.toFixed(2)}</div></div><div class="tile"><div class="k">time</div><div class="v small">${c.summary.elapsed_s} s</div></div></div>
  <p class="note small">Similarity is measured against the ion the chain started from, never against the previous round. The protected core of the original ion stays the reference throughout; a product already seen in this chain is never revisited. Click a step for its certificate.</p>
  <div class="chainrow"><div class="chainstep start"><div class="rnd"><span>start</span><span>${S.seed.heavy_atoms} atoms</span></div><div class="pic">${S.seed.svg}</div><div class="smi small mono" style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap">${esc(c.start.smiles)}</div></div>
  ${c.steps.map((s) => s.moved ? `<div class="chainstep" data-r="${s.round}"><div class="rnd"><span>round ${s.round}</span><span>${s.selected.similarity_to_start.toFixed(2)} to start</span></div><div class="pic">${s.selected.svg || ""}</div><div style="display:flex;justify-content:space-between;align-items:center">${opchip(s.selected.operator)}<span class="small">${s.candidates_certified} certified of ${s.candidates_examined}</span></div><div class="what small">${esc(s.selected.description)}</div></div>`
    : `<div class="chainstep" style="opacity:.6"><div class="rnd"><span>round ${s.round}</span></div><div class="empty" style="padding:20px">held: ${esc(s.reason)}</div></div>`).join("")}</div>
  <div class="actions"><span class="small">Download the chain</span><button class="btn sm" id="expChain">CSV</button></div></div>`;
  $$(".chainstep[data-r]", out).forEach((el) => el.addEventListener("click", () => { const s = c.steps.find((x) => x.round === Number(el.dataset.r)); openProduct({ ...s.selected, round: s.round, rule: c.rule }); }));
  $("#expChain").onclick = () => exportProducts(moved.map((s) => ({ ...s.selected, round: s.round, rule: c.rule })), "csv", "ionweaver_chain");
}

/* ---------------- batch ---------------- */
function renderBatch() {
  const root = $("#screen-batch");
  const sheet = S.batchSheet || null;
  const key = sheet ? sheet.sheet_id : "published";
  if (root.dataset.key === key) { if (S.batch && !$("#batchOut").innerHTML) paintBatch(); return; }
  root.dataset.key = key;
  const fams = sheet ? sheet.families.filter((f) => f.editable) : S.boot.families.filter((f) => f.editable);
  root.innerHTML = `<div class="kicker">Across families</div><h2 class="st">A quota of distinct ions per family, as one CSV</h2>
  <p class="lead">${sheet ? `The families of <b>${esc(sheet.filename)}</b>, edited against the published library plus the sheet's own components.` : `Every published family, or the ones you tick.`} Set the rules once and they hold for every family; no product appears twice, and none equals a parent.</p>
  <div class="two"><div class="card"><h3>Families <span class="small">· <a href="#" id="selAll">all</a> · <a href="#" id="selCat">cations</a> · <a href="#" id="selAn">anions</a> · <a href="#" id="selNone">none</a></span></h3>
    <div class="fam-select" id="famSel">${fams.map((f) => `<label><input type="checkbox" value="${esc(f.family_id)}" data-t="${f.ion_type}" checked> ${esc(f.family_label)} <span class="m">${f.editable} editable</span></label>`).join("")}</div>
    ${sheet ? `<div class="actions"><button class="btn sec sm" id="batchPub">Use the published families instead</button></div>` : (S.sheet ? `<div class="actions"><button class="btn sec sm" id="batchSheetBtn">Use my uploaded sheet instead</button></div>` : "")}</div>
  <div class="card">${rangeKnob("perFam", "Distinct ions per family", 1, 500, 25, 1, "")}${rangeKnob("ionsFam", "Parents used per family", 1, 40, 5, 1, "spread from the family's smallest editable ion to its largest; more parents, more variety")}
    <div class="field"><label>Operators</label><div class="chips" id="bops">${S.boot.operators.map((o) => `<span class="chip tog on" data-op="${o.key}">${o.key}</span>`).join("")}</div></div>
    ${constraintsPanel("consBatch")}
    <div class="actions"><button class="btn go" id="goBatch">Generate the set</button><span class="small" id="bEst"></span></div>
    <div id="batchProg" hidden><div class="progress"><i></i></div><div class="small pmsg"></div></div></div></div>
  <div id="batchOut"></div>`;
  wireRange("perFam", est); wireRange("ionsFam", est); wireConstraints("consBatch");
  const boxes = () => $$("#famSel input");
  const sel = (fn) => { boxes().forEach((b) => (b.checked = fn(b))); est(); };
  $("#selAll").onclick = (e) => { e.preventDefault(); sel(() => true); }; $("#selNone").onclick = (e) => { e.preventDefault(); sel(() => false); };
  $("#selCat").onclick = (e) => { e.preventDefault(); sel((b) => b.dataset.t === "cation"); }; $("#selAn").onclick = (e) => { e.preventDefault(); sel((b) => b.dataset.t === "anion"); };
  boxes().forEach((b) => b.addEventListener("change", est));
  $$("#bops .tog").forEach((c) => c.addEventListener("click", () => c.classList.toggle("on")));
  if ($("#batchPub")) $("#batchPub").onclick = () => { S.batchSheet = null; renderBatch(); };
  if ($("#batchSheetBtn")) $("#batchSheetBtn").onclick = () => { S.batchSheet = S.sheet; renderBatch(); };
  function est() { const n = boxes().filter((b) => b.checked).length; const per = Number($("#perFam").value); $("#bEst").textContent = `${n} famil${n === 1 ? "y" : "ies"} × ${per} = up to ${fmt(n * per)} ions`; }
  est();
  $("#goBatch").addEventListener("click", async () => {
    const ids = boxes().filter((b) => b.checked).map((b) => b.value); if (!ids.length) { toast("Tick at least one family"); return; }
    const ops = $$("#bops .tog.on").map((c) => c.dataset.op); if (!ops.length) { toast("Include at least one operator"); return; }
    const per = Number($("#perFam").value);
    const body = { family_ids: ids, per_family: per, ions_per_family: Number($("#ionsFam").value), operators: ops, constraints: readConstraints("consBatch"), svg: ids.length * per <= 400, sheet_id: sheet ? sheet.sheet_id : null };
    const prog = $("#batchProg"), btn = $("#goBatch"); prog.hidden = false; btn.disabled = true; $("#batchOut").innerHTML = "";
    const token = (S.batchToken = (S.batchToken || 0) + 1);
    try { const b = await job("/api/batch", body, (st) => progressBar(prog, st)); if (S.batchToken !== token) return; S.batch = b; paintBatch(); }
    catch (e) { if (S.batchToken === token) $("#batchOut").innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
    if (S.batchToken === token && $("#batchProg")) { $("#batchProg").hidden = true; $("#goBatch").disabled = false; }
  });
  if (S.batch) paintBatch();
}
function paintBatch() {
  const b = S.batch, out = $("#batchOut"), sm = b.summary;
  const ops = S.boot.operators.map((o) => o.key);
  out.innerHTML = `<div class="card"><h3>The set</h3>
  <div class="tiles"><div class="tile green"><div class="k">distinct certified ions</div><div class="v">${fmt(sm.delivered)}</div></div><div class="tile"><div class="k">families</div><div class="v">${sm.families}</div></div><div class="tile ink"><div class="k">new to the inventory</div><div class="v">${fmt(sm.new_to_inventory)}</div></div><div class="tile"><div class="k">time</div><div class="v small">${sm.elapsed_s} s</div></div></div>
  ${b.settings.constraints && Object.keys(b.settings.constraints).length ? `<p class="note"><b>Rules in force:</b> ${esc(constraintsText(b.settings.constraints))}.</p>` : ""}
  <div class="actions"><span class="small">Download the set</span><button class="btn sm" id="bCsv">CSV</button><button class="btn sec sm" id="bSdf">SDF</button><button class="btn sec sm" id="bSmi">SMILES</button></div>
  <div style="overflow-x:auto;margin-top:14px"><table><thead><tr><th>family</th><th></th><th class="num">asked</th><th class="num">delivered</th><th class="num">parents used</th><th class="num">new</th>${ops.map((o) => `<th class="num">${o}</th>`).join("")}</tr></thead>
  <tbody>${b.families.map((f, i) => `<tr class="rowin" style="animation-delay:${i * 20}ms"><td><b>${esc(f.family_label)}</b>${f.parents_breaking_rules ? `<br><span class="small" style="color:var(--amber)">${f.parents_breaking_rules} of ${f.ions_used} parents already break a body rule</span>` : ""}</td><td>${f.ion_type ? ionType(f.ion_type) : ""}</td><td class="num">${f.requested}</td><td class="num" style="color:${f.delivered < f.requested ? "var(--amber)" : "var(--green)"};font-weight:600">${f.delivered}</td><td class="num">${f.ions_used}${f.ions_available ? " of " + f.ions_available : ""}</td><td class="num">${f.new_to_inventory ?? ""}</td>${ops.map((o) => `<td class="num">${(f.operators || {})[o] ?? ""}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
  ${b.families.some((f) => f.delivered < f.requested) ? `<p class="note small">A family short of its quota had fewer distinct certified products available from its parents under these rules. A parent whose own body already breaks a rule (a fluorine on a BF₄ derivative under "carbon only") can only comply through an edit that removes the offending atoms. Raise the number of parents, loosen the rule for that family, or accept the smaller set.</p>` : ""}</div>
  ${b.products[0] && b.products[0].svg ? `<div class="card"><h3>A look at them</h3><div class="pgrid" id="bgrid">${b.products.map((p, i) => productCard(p, i)).join("")}</div></div>` : `<p class="note small">Previews are drawn for sets of up to 400 ions; this set is in the CSV.</p>`}`;
  const name = `ionweaver_${b.settings.workspace === "published" ? "families" : "sheet"}_${sm.delivered}`;
  $("#bCsv").onclick = () => exportProducts(b.products, "csv", name); $("#bSdf").onclick = () => exportProducts(b.products, "sdf", name); $("#bSmi").onclick = () => exportProducts(b.products, "smi", name);
  $$("#bgrid .pcard").forEach((c) => c.addEventListener("click", () => openProduct(b.products[Number(c.dataset.i)])));
}

/* ---------------- fragment library ---------------- */
async function renderLibrary() {
  const root = $("#screen-library");
  if (root.dataset.ready) return;
  root.innerHTML = `<div class="kicker">Component library</div><h2 class="st">The side groups the curated ions display</h2><div class="progress busy"><i style="width:30%"></i></div>`;
  let d;
  try { d = await api("/api/library"); } catch (e) { root.innerHTML = `<div class="note bad">${esc(e.message)}</div>`; return; }
  root.dataset.ready = "1";
  const sm = d.summary;
  root.innerHTML = `<div class="kicker">Component library</div><h2 class="st">The side groups the curated ions display</h2>
  <p class="lead">Every side group of every editable ion, cut at the bond that joins it to the core and carried whole with that bond's order; the dummy atom marks the attachment. Graft attaches a component at a free valence, Swap exchanges a whole side group for one, and Cross-over splices on one of the ${d.fragments} tails cut from these. A one-atom component is left to Insert and Substitute.</p>
  <div class="tiles"><div class="tile ink"><div class="k">components</div><div class="v">${sm.library_components}</div></div><div class="tile"><div class="k">instances across the ions</div><div class="v">${fmt(sm.component_instances)}</div></div><div class="tile"><div class="k">eligible for Graft</div><div class="v">${sm.eligible.graft}</div></div><div class="tile"><div class="k">eligible for Swap</div><div class="v">${sm.eligible.swap}</div></div><div class="tile"><div class="k">tails for Cross-over</div><div class="v">${d.fragments}</div></div><div class="tile ink"><div class="k">largest, in atoms</div><div class="v">${sm.largest_component_heavy_atoms}</div></div></div>
  <div class="filter-row" style="margin-top:14px"><span class="chip tog on" data-f="all">All</span><span class="chip tog" data-f="graft">Graft-eligible</span><span class="chip tog" data-f="ring">With a ring</span><span class="chip tog" data-f="hetero">With a heteroatom</span><input type="text" id="libSearch" placeholder="search SMILES or a family…"><select id="libSort" style="width:auto"><option value="occ">most common first</option><option value="size_asc">smallest first</option><option value="size_desc">largest first</option></select><span class="small" id="libCount"></span></div>
  <div class="grid" id="libGrid"></div>`;
  let f = "all", q = "", sort = "occ";
  const paint = () => {
    let rows = d.components.filter((c) => (f === "all" || (f === "graft" && c.graft) || (f === "ring" && c.rings > 0) || (f === "hetero" && /[NOSPFBClBrI]/.test(c.key.replace(/\*/g, "")))) && (c.key.toLowerCase().includes(q) || c.family_labels.join(" ").toLowerCase().includes(q)));
    rows = rows.sort((a, b) => sort === "occ" ? b.occurrences - a.occurrences : sort === "size_asc" ? a.heavy_atoms - b.heavy_atoms : b.heavy_atoms - a.heavy_atoms);
    $("#libCount").textContent = `${rows.length} shown`;
    $("#libGrid").innerHTML = rows.map((c, i) => `<div class="gcard" style="cursor:default;animation-delay:${Math.min(i * 10, 400)}ms"><div class="pic" style="min-height:120px">${c.svg}</div>
      <div class="ttl" style="font-size:14px"><span class="mono">${esc(c.key)}</span><span class="small">${c.library_id}</span></div>
      <div class="meta"><span>${c.heavy_atoms} atom${c.heavy_atoms === 1 ? "" : "s"}</span><span>${c.occurrences} occurrence${c.occurrences === 1 ? "" : "s"}</span><span>${c.ions} ion${c.ions === 1 ? "" : "s"}</span><span>${c.families} famil${c.families === 1 ? "y" : "ies"}</span></div>
      <div class="chips" style="margin-top:4px">${c.graft ? '<span class="chip" style="padding:2px 9px;font-size:11.5px">Graft</span>' : ""}${c.swap ? '<span class="chip" style="padding:2px 9px;font-size:11.5px">Swap</span>' : ""}${c.crossover ? '<span class="chip" style="padding:2px 9px;font-size:11.5px">Cross-over</span>' : ""}${c.prune ? '<span class="chip" style="padding:2px 9px;font-size:11.5px">Prune</span>' : ""}${!c.graft && !c.swap ? `<span class="chip grey" style="padding:2px 9px;font-size:11.5px" title="${esc(c.reason)}">Insert or Substitute only</span>` : ""}</div>
      <div class="small" style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="${esc(c.family_labels.join(", "))}">${esc(c.family_labels.slice(0, 4).join(", "))}${c.family_labels.length > 4 ? "…" : ""}</div></div>`).join("");
  };
  $$(".filter-row .tog", root).forEach((c) => c.addEventListener("click", () => { $$(".filter-row .tog", root).forEach((x) => x.classList.remove("on")); c.classList.add("on"); f = c.dataset.f; paint(); }));
  $("#libSearch").addEventListener("input", (e) => { q = e.target.value.trim().toLowerCase(); paint(); });
  $("#libSort").addEventListener("change", (e) => { sort = e.target.value; paint(); });
  paint();
}

/* ---------------- about ---------------- */
function renderAbout() {
  const b = S.boot, root = $("#screen-about"), e = b.engine, h = b.headline;
  root.innerHTML = `<div class="kicker">Method &amp; numbers</div><h2 class="st">What IonWeaver does, and what it does not claim</h2>
  <div class="card"><h3>The idea</h3><p class="note" style="margin-top:0">A family of ions is defined by a charged group. IonWeaver finds the largest substructure that every member of a family shares, carrying its charged atom, and calls it the protected core. An edit is <b>defined on the parent graph</b> (what it removes, what it puts back, where) and decided there against a valence ledger and the sets of atom and bond environments the curated inventory displays; only then is a product built. Every product is <b>certified</b> against five requirements before it is offered. Across the submitted study, certification refused ${h.products_refused_by_certification} of ${fmt(h.products_certified_at_the_reported_allowance)} products: the construction rules held, which is the claim.</p></div>
  <div class="two"><div class="card"><h3>The eight operators</h3><table><tbody>${b.operators.map((o) => `<tr><td><b>${o.key}</b><br><span class="small">${o.level}-level</span></td><td>${esc(o.text)}</td></tr>`).join("")}</tbody></table></div>
  <div class="card"><h3>The five requirements, tested on every product</h3><div class="checks">${b.requirements.map((r, i) => `<div class="check"><i>${i + 1}</i><div><b>${esc(r.short)}.</b> ${esc(r.text)}</div></div>`).join("")}</div>
  <p class="note small">Requirements 3 and 4 read labels applied to the core before the edit, so no substructure search is run on a product: a product that moved the core cannot certify.</p>
  <h3 style="margin-top:16px">Conditions on a definition, before anything is built</h3><p class="note small" style="margin-top:0">The valence ledger (every touched atom lands on a legal valence), the atom-environment condition (every touched atom lands in an environment the inventory displays), the bond-environment condition (every new or changed bond is a pair the inventory displays), and three restrictions on constructed atoms (neutral and non-aromatic; no terminal halogen; not above the lowest valence that admits hydrogen).</p></div></div>
  <div class="card"><h3>The submitted run</h3><div class="tiles">
    <div class="tile"><div class="k">source entries read</div><div class="v">${fmt(h.source_entries_read)}</div></div><div class="tile"><div class="k">ions retained</div><div class="v">${fmt(h.retained_ions)}</div></div><div class="tile"><div class="k">families</div><div class="v">${h.families_retained}</div></div><div class="tile"><div class="k">cores accepted</div><div class="v">${h.cores_accepted}</div></div><div class="tile"><div class="k">ions a core reached</div><div class="v">${fmt(h.ions_a_core_reached)}</div></div>
    <div class="tile purple"><div class="k">library components</div><div class="v">${h.library_components}</div></div><div class="tile purple"><div class="k">terminal fragments</div><div class="v">${h.terminal_fragments_at_two_or_more_heavy_atoms}</div></div><div class="tile amber"><div class="k">definitions offered</div><div class="v">${fmt(h.definitions_offered)}</div></div><div class="tile amber"><div class="k">definitions accepted</div><div class="v">${fmt(h.definitions_accepted)}</div></div>
    <div class="tile green"><div class="k">examinations, allowance 100</div><div class="v">${fmt(h.examinations)}</div></div><div class="tile green"><div class="k">certified, allowance 50</div><div class="v">${fmt(h.products_certified_at_the_reported_allowance)}</div></div><div class="tile green"><div class="k">distinct products</div><div class="v">${fmt(h.distinct_products_at_the_reported_allowance)}</div></div><div class="tile ink"><div class="k">signature collisions</div><div class="v">${h.signature_collisions}</div></div><div class="tile ink"><div class="k">perception mismatches</div><div class="v">${h.perception_check_mismatches}</div></div>
    <div class="tile green"><div class="k">validity</div><div class="v">${pct(h.validity_percent)}</div></div><div class="tile"><div class="k">uniqueness</div><div class="v">${pct(h.uniqueness_percent)}</div></div><div class="tile"><div class="k">novelty</div><div class="v">${pct(h.novelty_percent)}</div></div><div class="tile"><div class="k">internal diversity</div><div class="v">${Number(h.internal_diversity).toFixed(3)}</div></div><div class="tile"><div class="k">median similarity, product to parent</div><div class="v">${Number(h.median_similarity_between_a_product_and_the_ion_it_came_from).toFixed(3)}</div></div></div>
  <div class="two" style="margin-top:16px"><div><h3>Allowance sensitivity</h3><table><thead><tr><th>allowance</th><th class="num">examinations</th><th class="num">certified</th><th class="num">distinct</th><th class="num">ions reached</th></tr></thead><tbody>${b.sensitivity.map((r) => `<tr><td>${r.allowance}</td><td class="num">${fmt(r.examinations)}</td><td class="num">${fmt(r.certified)}</td><td class="num">${fmt(r.distinct_products)}</td><td class="num">${fmt(r.ions_reached)}</td></tr>`).join("")}</tbody></table></div>
  <div><h3>Matched string controls (one elementary string edit per graph candidate)</h3><table><thead><tr><th>method</th><th class="num">attempts</th><th class="num">accepted</th><th class="num">share</th></tr></thead><tbody>${b.string_controls.map((r) => `<tr><td>${r.method.toUpperCase()}</td><td class="num">${fmt(r.attempts)}</td><td class="num">${fmt(r.accepted)}</td><td class="num">${pct(100 * r.accepted / r.attempts)}</td></tr>`).join("")}<tr><td><b>IonWeaver</b></td><td class="num">${fmt(h.examinations_at_the_reported_allowance)}</td><td class="num">${fmt(h.products_certified_at_the_reported_allowance)}</td><td class="num">${pct(100 * h.products_certified_at_the_reported_allowance / h.examinations_at_the_reported_allowance)}</td></tr></tbody></table></div></div></div>
  <div class="card"><h3>The engine running here</h3><p class="note" style="margin-top:0">${e.unmodified ? `<span class="chip green">unmodified</span> every module matches the submitted notebook (<code>${esc(e.source_notebook)}</code>, vendored ${esc(e.vendored_utc)}).` : `<span class="chip red">modified</span> at least one module differs from the notebook.`} RDKit ${esc(e.pins.RDKIT_VERSION)} pinned by the method; run seed ${esc(e.pins.SEED)}. The library rebuilt at start-up reproduced the archive's curation, registries, editable set and component library: ${Object.entries(b.library_checks).filter(([k]) => k !== "all").every(([, v]) => v) ? "every check held" : "a check failed"} (${b.load_seconds} s).</p>
  <table><thead><tr><th>module</th><th class="num">lines</th><th>sha-256</th><th></th></tr></thead><tbody>${Object.entries(e.modules).map(([k, m]) => `<tr><td class="mono">iw/${k}.py</td><td class="num">${m.lines}</td><td class="mono small">${m.sha256.slice(0, 24)}…</td><td>${m.unmodified ? '<span class="chip green" style="padding:2px 9px;font-size:12px">matches</span>' : '<span class="chip red" style="padding:2px 9px;font-size:12px">differs</span>'}</td></tr>`).join("")}</tbody></table></div>
  <div class="card"><h3>What it does not claim</h3><p class="note" style="margin-top:0">Each output is one cation or one anion, not an ionic liquid: pairing with a counter-ion, property prediction and synthesis assessment all come afterwards. An ion absent from the reference inventory is not thereby unknown in the wider literature. Stereochemistry is not preserved or evaluated. Core discovery is automatic, but the initial family grouping is not: IonWeaver needs labelled examples and does not classify an unlabelled collection into families. The descriptors shown are RDKit's atom-contribution estimates on a charged ion and are offered as filters, not as predictions.</p></div>`;
}

/* ---------------- boot ---------------- */
const RENDER = { start: renderStart, families: renderFamilies, own: renderOwn, sheet: renderSheet, bench: renderBench, batch: renderBatch, library: renderLibrary, about: renderAbout, studio: () => renderStudio() };
async function boot() {
  // a free host puts the engine to sleep between visits and takes up to a minute to wake
  const started = Date.now();
  while (!S.boot) {
    try { S.boot = await api("/api/bootstrap"); }
    catch (err) {
      if (Date.now() - started > 120000) { $("#engineDot").className = "dot bad"; $("#engineText").textContent = "engine unreachable"; $("#screen-start").innerHTML = `<div class="note bad"><b>The engine server could not be reached</b> (${esc(API || "same origin")}). ${esc(err.message)}<br>Reload in a minute; if it persists, the server is down.</div>`; return; }
      $("#engineText").textContent = "waking the engine…";
      $("#screen-start").innerHTML = `<div class="hero"><div class="kicker">One moment</div><h1 class="t">Waking the engine…</h1><p class="lead">The engine server sleeps between visits and takes up to a minute to start. This page will continue on its own.</p><div class="progress busy" style="max-width:420px"><i style="width:30%"></i></div></div>`;
      await new Promise((r) => setTimeout(r, 4000));
    }
  }
  try {
    const e = S.boot.engine, c = S.boot.counts;
    $("#engineDot").className = "dot " + (e.unmodified && S.boot.library_checks.all ? "ok" : "bad");
    $("#engineText").textContent = e.unmodified ? "engine · unmodified" : "engine · modified";
    $("#engineBadge").onclick = () => show("about");
    $("#navNote").innerHTML = `${fmt(c.ions)} ions · ${c.families} families · ${c.cores} cores<br>RDKit ${esc(e.pins.RDKIT_VERSION)} · served ${API ? "by " + esc(API.replace(/^https?:\/\//, "")) : "locally"}`;
    const hash = location.hash.replace("#", "");
    show(RENDER[hash] ? hash : "start");
  } catch (err) {
    $("#engineDot").className = "dot bad"; $("#engineText").textContent = "page error";
    $("#screen-start").innerHTML = `<div class="note bad"><b>The page could not start.</b> ${esc(err.message)}</div>`;
  }
}
boot();
