/* IonWeaver · the studio. Build an ion one chosen edit at a time. Every verdict here comes from
   the engine: a move is either among the definitions it offers (then it is built and certified)
   or the engine names the rule that keeps it out. */
"use strict";

const ORDER_GLYPH = { SINGLE: "–", DOUBLE: "=", TRIPLE: "≡" };
const ORDER_WORD = { SINGLE: "single", DOUBLE: "double", TRIPLE: "triple" };

async function openStudio(start) {
  if (!start) { show("studio"); return; }
  show("studio");
  const root = $("#screen-studio");
  root.innerHTML = `<div class="kicker">Studio</div><h2 class="st">Opening…</h2><div class="progress busy"><i style="width:30%"></i></div><p class="small">deciding every move the eight operators offer on this ion</p>`;
  try {
    const cat = await api("/api/studio/open", start);
    const keep = S.studio || {};
    S.studio = { trail: [{ state: cat.state, record: null, cat }], site: null, verdict: null, pieces: keep.pieces || {}, libData: keep.libData || null, libPromise: keep.libPromise || null, rule: keep.rule || "random" };
    renderStudio();
  } catch (e) { root.innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
}
const cur = () => S.studio.trail[S.studio.trail.length - 1];

function renderStudio() {
  const root = $("#screen-studio");
  if (!S.studio) {
    root.innerHTML = `<div class="kicker">Studio</div><h2 class="st">Build an ion by hand, one certified edit at a time</h2>
    <p class="lead">Point at an atom or a bond, choose a component, an element or a cut, and the engine either makes the move and certifies the result or names the rule that stops it. Every step becomes the new parent; the original core stays protected throughout.</p>
    <div class="card"><h3>Start from</h3><div class="actions" style="margin-top:0">
      ${S.start ? `<button class="btn" id="stFromBench">The workbench ion (${esc(startName())})</button>` : ""}
      <button class="btn sec" data-go="families">A published family</button><button class="btn sec" data-go="own">Your own ion</button><button class="btn sec" data-go="sheet">Your own family</button>
      <button class="btn sec" id="stExample">Try it on 1-butyl-3-methylimidazolium</button></div></div>`;
    $$("[data-go]", root).forEach((b) => b.addEventListener("click", () => show(b.dataset.go)));
    if ($("#stFromBench")) $("#stFromBench").onclick = () => openStudio(S.start);
    $("#stExample").onclick = () => openStudio({ kind: "published", inventory_id: 42 });
    return;
  }
  const st = S.studio, node = cur(), cat = node.cat, state = node.state, m = cat.moves;
  const atomMoves = (i) => (m.graft[i] || []).length + (m.append[i] || []).length + (m.substitute[i] || []).length + (m.delete.includes(i) ? 1 : 0) + cat.components.filter((c) => c.anchor === i).reduce((n, c) => n + (m.swap[c.idx] || []).length, 0);
  const bondKey = (a, b) => (a < b ? `${a}-${b}` : `${b}-${a}`);
  const bondMoves = (a, b) => (m.insert[bondKey(a, b)] || []).length + (m.rebond[bondKey(a, b)] || []).length + (m.prune.some((p) => (p.anchor === a && p.root === b) || (p.anchor === b && p.root === a)) ? 1 : 0) + Object.values(m.crossover).reduce((n, cuts) => n + Object.values(cuts).filter((c) => (c.anchor === a && c.root === b) || (c.anchor === b && c.root === a)).reduce((k, c) => k + c.fragments.length, 0), 0);
  const last = node.record;
  const siteOn = (t, ...xs) => st.site && st.site.type === t && xs.every((x, k) => st.site.v[k] === x);
  root.innerHTML = `<div class="kicker">Studio · ${esc(state.workspace === "published" ? "published inventory" : "published inventory + your sheet")}</div>
  <div class="results-head"><h2 class="st" style="margin:0">${esc(state.family_label || state.family_id)} ${ionType(state.family_id.startsWith("cation") ? "cation" : "anion")} <span class="small">· core ${state.grammar_id} · step ${state.step}</span></h2>
  <div class="actions" style="margin:0"><select id="pickRule" style="width:auto">${S.boot.rules.map((r) => `<option value="${r.key}" ${st.rule === r.key ? "selected" : ""}>${esc(r.label)}</option>`).join("")}</select><button class="btn go sm" id="pickBtn">Let the engine choose</button><button class="btn sec sm" id="undoBtn" ${st.trail.length > 1 ? "" : "disabled"} title="${st.trail.length > 1 ? "" : "nothing to undo yet"}">Undo</button><button class="btn sec sm" id="toBench">Open in the workbench</button><button class="btn sec sm" id="trailCsv" ${st.trail.length > 1 ? "" : "disabled"} title="${st.trail.length > 1 ? "" : "no step yet"}">Download the trail (CSV)</button></div></div>
  <div class="stage" style="margin-top:12px">
    <div class="card"><h3>The ion now</h3><div class="small" style="margin:-6px 0 8px">Body atoms are numbered; the core is labelled Q.</div>
      <div class="sitepic" id="sitepic"><span class="stepno">step ${state.step}</span>${cat.svg}</div>
      <dl class="facts"><dt>SMILES</dt><dd class="mono">${esc(state.smiles)}</dd><dt>size</dt><dd>${cat.heavy_atoms} atoms (H not counted) · ${fmt(cat.total)} moves on offer</dd>${last ? `<dt>last move</dt><dd>${opchip(last.operator)} ${esc(last.description)}</dd><dt>similarity to the start</dt><dd>${last.similarity_to_origin !== null ? last.similarity_to_origin.toFixed(2) : "–"} · ${last.new_to_inventory ? "new to the inventory" : "in the inventory"}</dd>` : ""}</dl>
      <h3 style="margin-top:14px">Pick a site</h3>
      <div class="sites">${cat.atoms.map((a) => `<span class="site ${a.core ? "core" : ""} ${siteOn("atom", a.idx) ? "on" : ""} ${atomMoves(a.idx) ? "" : "none"}" data-atom="${a.idx}" title="${a.core ? "core atom Q" + a.q : "body atom " + a.idx}${a.fixed && !a.core ? " · held in a ring with the core" : ""}">${a.core ? "Q" + a.q : a.idx} <small>${a.symbol}${a.hs ? "H" + (a.hs > 1 ? a.hs : "") : ""} · ${atomMoves(a.idx)}</small></span>`).join("")}</div>
      <div class="sites">${cat.bonds.filter((b) => b.editable).map((b) => `<span class="site ${siteOn("bond", b.a, b.b) ? "on" : ""} ${bondMoves(b.a, b.b) ? "" : "none"}" data-bond="${b.a},${b.b}" title="bond ${b.a}${ORDER_GLYPH[b.order] || "-"}${b.b}">${lbl(cat, b.a)}${ORDER_GLYPH[b.order] || "–"}${lbl(cat, b.b)} <small>${bondMoves(b.a, b.b)}</small></span>`).join("")}</div>
      <p class="small" style="margin-top:8px">The number is how many moves the engine offers there. A site with none is still shown so you can ask why.</p></div>
    <div class="card" id="movePanel">${st.site ? "" : `<div class="empty">Pick an atom or a bond on the left.<br><span class="small">Or let the engine choose the next move.</span></div>`}</div>
  </div>
  <div class="card"><h3>Trail</h3><div class="small" style="margin:-6px 0 8px">Click a step to return to it.</div><div class="trail" id="trail">${st.trail.map((n, i) => n.record ? `<div class="chainstep ${i === st.trail.length - 1 ? "cur" : ""}" data-i="${i}"><div class="rnd"><span>step ${i}</span><span>${n.record.similarity_to_origin !== null ? n.record.similarity_to_origin.toFixed(2) + " to start" : ""}</span></div><div class="pic">${n.record.svg || ""}</div><div style="display:flex;justify-content:space-between;align-items:center">${opchip(n.record.operator)}<span class="small">${n.record.descriptors.formula}</span></div><div class="what">${esc(n.record.description)}</div></div>`
    : `<div class="chainstep start ${i === st.trail.length - 1 ? "cur" : ""}" data-i="0"><div class="rnd"><span>start</span><span>${n.cat.heavy_atoms} atoms</span></div><div class="pic">${n.cat.svg}</div><div class="smi small mono" style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap">${esc(n.state.smiles)}</div></div>`).join("")}</div></div>`;
  $$(".sites [data-atom]", root).forEach((el) => el.addEventListener("click", () => selectSite({ type: "atom", v: [Number(el.dataset.atom)] })));
  $$(".sites [data-bond]", root).forEach((el) => el.addEventListener("click", () => { const [a, b] = el.dataset.bond.split(",").map(Number); selectSite({ type: "bond", v: [a, b] }); }));
  $("#pickRule").addEventListener("change", (e) => (st.rule = e.target.value));
  $("#pickBtn").addEventListener("click", enginePick);
  $("#undoBtn").addEventListener("click", () => goToStep(st.trail.length - 2));
  $("#toBench").addEventListener("click", () => openBench({ kind: "state", smiles: state.smiles, core_atoms: state.core_atoms, family_id: state.family_id, family_label: state.family_label, grammar_id: state.grammar_id, pattern: state.core_pattern, workspace: state.workspace, inventory_id: state.step ? -1 : state.inventory_id }));
  $("#trailCsv").addEventListener("click", () => exportProducts(st.trail.filter((n) => n.record).map((n, i) => ({ ...n.record, round: i + 1, rule: "studio" })), "csv", "ionweaver_studio_trail"));
  $$("#trail .chainstep", root).forEach((el) => el.addEventListener("click", () => goToStep(Number(el.dataset.i))));
  if (st.site) paintMovePanel();
}
function lbl(cat, i) { const a = cat.atoms[i]; return a.core ? "Q" + a.q : String(i); }
function goToStep(i) { const st = S.studio; if (i < 0 || i >= st.trail.length) return; st.trail = st.trail.slice(0, i + 1); st.site = null; st.verdict = null; renderStudio(); }

async function selectSite(site) {
  const st = S.studio; st.site = site; st.verdict = null;
  $$(".sites .site").forEach((el) => el.classList.remove("on"));
  const sel = site.type === "atom" ? $(`.sites [data-atom="${site.v[0]}"]`) : $(`.sites [data-bond="${site.v[0]},${site.v[1]}"]`);
  if (sel) sel.classList.add("on");
  const focus = api("/api/studio/focus", { state: cur().state, focus: site.v }).then((d) => { const pic = $("#sitepic"); if (pic && st.site === site) pic.innerHTML = `<span class="stepno">step ${cur().state.step}</span>` + d.svg; }).catch(() => {});
  await paintMovePanel();
  await focus;
}

async function pieces(keys) {
  const st = S.studio;
  if (!st.libPromise) st.libPromise = api("/api/library").then((d) => { d.components.forEach((c) => (st.pieces[c.key] = c.svg)); return d; }).catch(() => ({ components: [] }));
  st.libData = await st.libPromise;
  const still = keys.filter((k) => !(k in st.pieces));
  if (still.length) { try { Object.assign(st.pieces, await api("/api/piece", { keys: still.slice(0, 600) })); } catch (_) {} }
}

function elemChips(items, allowed, mk, withOrder) {
  // items: [[el, order]] or [el]; allowed: Set of "el|order" or "el"
  return `<div class="elems">${items.map((it) => { const el = withOrder ? it[0] : it, order = withOrder ? it[1] : null, key = withOrder ? `${el}|${order}` : el, ok = allowed.has(key); return `<button type="button" class="elem ${ok ? "" : "no"}" data-move='${esc(JSON.stringify(mk(el, order)))}' title="${ok ? "offered" : "not offered here: click to ask why"}">${el}${withOrder ? `<small>${ORDER_WORD[order]}</small>` : ""}</button>`; }).join("")}</div>`;
}
function fragGrid(id, entries, allowedSet, mk, hint) {
  return `<div class="field" style="margin:6px 0"><input type="text" class="fragsearch" data-for="${id}" placeholder="search ${hint}…"></div>
  <div class="fraggrid" id="${id}">${entries.map((e) => `<div class="frag ${allowedSet.has(e.key) ? "" : "no"}" data-key="${esc(e.key)}" data-move='${esc(JSON.stringify(mk(e.key)))}' title="${allowedSet.has(e.key) ? "offered" : "not offered here: click to ask why"}"><div class="pic">${S.studio.pieces[e.key] || ""}</div><div class="k">${esc(e.key)}</div><div class="n">${e.heavy_atoms} atoms${e.occurrences ? " · " + e.occurrences + "×" : ""}</div></div>`).join("")}</div>`;
}

async function paintMovePanel() {
  const st = S.studio, node = cur(), cat = node.cat, m = cat.moves, site = st.site, panel = $("#movePanel");
  if (!site) return;
  const atomsIn = (idx) => { const a = cat.atoms[idx]; return `${a.core ? "core atom Q" + a.q : "atom " + idx} (${a.symbol}${a.hs ? "H" + (a.hs > 1 ? a.hs : "") : ""}${a.ring ? ", in a ring" : ""}${a.aromatic ? ", aromatic" : ""})`; };
  const blocks = [];
  const graftLib = cat.library.filter((e) => e.graft).map((e) => ({ ...e, occurrences: (st.libData ? (st.libData.components.find((c) => c.key === e.key) || {}).occurrences : 0) }));
  if (site.type === "atom") {
    const i = site.v[0], a = cat.atoms[i];
    panel.innerHTML = `<h3>At ${atomsIn(i)}</h3><div class="progress busy"><i style="width:30%"></i></div>`;
    await pieces(graftLib.map((e) => e.key));
    const allowedG = new Set(m.graft[String(i)] || []);
    const libSorted = graftLib.map((e) => ({ ...e, occurrences: (st.libData.components.find((c) => c.key === e.key) || {}).occurrences || 0 })).sort((x, y) => (allowedG.has(y.key) - allowedG.has(x.key)) || (y.occurrences - x.occurrences) || (x.heavy_atoms - y.heavy_atoms));
    blocks.push(`<div class="moveblock"><h4>Attach a component <span class="small">Graft · ${allowedG.size} of ${graftLib.length} offered here</span></h4>${fragGrid("gGrid", libSorted, allowedG, (k) => ({ operator: "Graft", atom: i, structure: k }), "a component")}</div>`);
    const allowedA = new Set((m.append[String(i)] || []).map(([e, o]) => `${e}|${o}`));
    const items = []; cat.elements.forEach((e) => cat.orders.forEach((o) => items.push([e, o])));
    blocks.push(`<div class="moveblock"><h4>Append one atom <span class="small">Insert · ${allowedA.size} offered</span></h4>${elemChips(items, allowedA, (e, o) => ({ operator: "Insert", mode: "appended", atom: i, element: e, order: o }), true)}</div>`);
    const allowedS = new Set(m.substitute[String(i)] || []);
    blocks.push(`<div class="moveblock"><h4>Change the element <span class="small">Substitute · ${allowedS.size} offered</span></h4>${elemChips(cat.elements.filter((e) => e !== a.symbol), allowedS, (e) => ({ operator: "Substitute", atom: i, element: e }), false)}</div>`);
    const canDel = m.delete.includes(i);
    blocks.push(`<div class="moveblock"><h4>Remove this atom <span class="small">Delete</span></h4><button type="button" class="elem ${canDel ? "" : "no"}" data-move='${esc(JSON.stringify({ operator: "Delete", atom: i }))}'>remove<small>${a.degree === 2 ? "and rejoin" : a.degree === 1 ? "terminal" : "degree " + a.degree}</small></button></div>`);
    const comps = cat.components.filter((c) => c.anchor === i);
    for (const c of comps) {
      const swapLib = cat.library.filter((e) => e.swap && e.key !== c.key);
      await pieces(swapLib.map((e) => e.key));
      const allowedW = new Set(m.swap[String(c.idx)] || []);
      const sorted = swapLib.map((e) => ({ ...e, occurrences: (st.libData.components.find((x) => x.key === e.key) || {}).occurrences || 0 })).sort((x, y) => (allowedW.has(y.key) - allowedW.has(x.key)) || (y.occurrences - x.occurrences));
      blocks.push(`<div class="moveblock"><h4>Exchange the side group <code>${esc(c.key)}</code> (atoms ${c.atoms.join(", ")}) <span class="small">Swap · ${allowedW.size} offered${c.heavy_atoms < 2 ? " · a one-atom group is a Substitute" : ""}</span></h4>${fragGrid("sGrid" + c.idx, sorted, allowedW, (k) => ({ operator: "Swap", anchor: i, root: c.root, structure: k }), "a component")}</div>`);
    }
  } else {
    const [a, b] = site.v, key = a < b ? `${a}-${b}` : `${b}-${a}`, bond = cat.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    panel.innerHTML = `<h3>At the ${ORDER_WORD[bond.order] || ""} bond ${lbl(cat, a)}${ORDER_GLYPH[bond.order] || "–"}${lbl(cat, b)}</h3>`;
    const allowedI = new Set(m.insert[key] || []);
    blocks.push(`<div class="moveblock"><h4>Insert one atom into the bond <span class="small">Insert · ${allowedI.size} offered</span></h4>${elemChips(cat.elements, allowedI, (e) => ({ operator: "Insert", mode: "bridged", bond: [a, b], element: e }), false)}</div>`);
    const allowedR = new Set(m.rebond[key] || []);
    blocks.push(`<div class="moveblock"><h4>Change the bond order <span class="small">Rebond</span></h4><div class="elems">${cat.orders.filter((o) => o !== bond.order).map((o) => `<button type="button" class="elem ${allowedR.has(o) ? "" : "no"}" data-move='${esc(JSON.stringify({ operator: "Rebond", bond: [a, b], order: o }))}'>${ORDER_GLYPH[o]}<small>${ORDER_WORD[o]}</small></button>`).join("")}</div></div>`);
    const prunes = m.prune.filter((p) => (p.anchor === a && p.root === b) || (p.anchor === b && p.root === a));
    const cuts = [];
    Object.values(m.crossover).forEach((slot) => Object.values(slot).forEach((c) => { if ((c.anchor === a && c.root === b) || (c.anchor === b && c.root === a)) cuts.push(c); }));
    const compHere = cat.components.find((c) => c.atoms.includes(a) && c.atoms.includes(b));
    if (compHere) {
      blocks.push(`<div class="moveblock"><h4>Cut the piece beyond this bond <span class="small">Prune</span></h4>${prunes.length ? prunes.map((p) => `<button type="button" class="elem" data-move='${esc(JSON.stringify({ operator: "Prune", anchor: p.anchor, root: p.root }))}'>cut<small>${p.atoms} atoms away</small></button>`).join("") : `<button type="button" class="elem no" data-move='${esc(JSON.stringify({ operator: "Prune", anchor: a, root: b }))}'>cut<small>not offered</small></button>`}</div>`);
      const anchor = cuts.length ? cuts[0].anchor : a, rootAt = cuts.length ? cuts[0].root : b;
      const allowedX = new Set(cuts.flatMap((c) => c.fragments));
      const tails = cat.fragments.filter((f) => f.order === bond.order || allowedX.has(f.key));
      await pieces(tails.map((f) => f.key));
      const sortedT = tails.sort((x, y) => (allowedX.has(y.key) - allowedX.has(x.key)) || (x.heavy_atoms - y.heavy_atoms));
      blocks.push(`<div class="moveblock"><h4>Splice a different tail beyond this bond <span class="small">Cross-over · ${allowedX.size} of ${tails.length} tails offered</span></h4>${cuts.length ? `<p class="small" style="margin:0 0 6px">${cuts[0].removed_atoms} atoms beyond the bond would be replaced.</p>` : `<p class="small" style="margin:0 0 6px">No terminal piece lies beyond this bond; click a tail to see why.</p>`}${fragGrid("xGrid", sortedT, allowedX, (k) => ({ operator: "Cross-over", anchor, root: rootAt, structure: k }), "a tail")}</div>`);
    }
  }
  if (S.studio !== st || st.site !== site || cur() !== node) return;   // the user moved on while we were drawing
  panel.innerHTML = `<h3>${site.type === "atom" ? "At " + atomsIn(site.v[0]) : panel.querySelector("h3").innerHTML}</h3><div id="verdict"></div>${blocks.join("")}`;
  $$("[data-move]", panel).forEach((el) => el.addEventListener("click", () => tryMove(JSON.parse(el.dataset.move))));
  $$(".fragsearch", panel).forEach((inp) => inp.addEventListener("input", () => { const q = inp.value.trim().toLowerCase(); $$(`#${inp.dataset.for} .frag`).forEach((f) => (f.hidden = q && !f.dataset.key.toLowerCase().includes(q))); }));
  if (st.verdict) paintVerdict(st.verdict);
}
function paintVerdict(v) {
  const el = $("#verdict"); if (!el) return;
  el.innerHTML = v.ok ? `<div class="verdict yes"><b>Done.</b> ${esc(v.description || "")}</div>`
    : `<div class="verdict no"><b>Not here.</b> ${esc(v.reason || "")}${v.condition ? ` <span class="small">(${esc(v.condition)}${v.at ? ", " + esc(v.at) : ""}${v.environment ? ": " + esc(v.environment) : ""})</span>` : ""}</div>`;
}

async function tryMove(move) {
  const st = S.studio, node = cur();
  const el = $("#verdict"); if (el) el.innerHTML = `<div class="verdict yes" style="background:#EEF3FA;color:var(--ku2);border-color:var(--ku2)">Asking the engine…</div>`;
  try {
    const r = await api("/api/studio/try", { state: node.state, move });
    if (S.studio !== st || cur() !== node) return;          // the trail changed while the engine worked
    if (!r.ok) { st.verdict = r; paintVerdict(r); return; }
    st.trail.push({ state: r.state, record: r.record, cat: r.catalogue }); st.site = null; st.verdict = null;
    renderStudio();
  } catch (e) { st.verdict = { ok: false, reason: e.message }; paintVerdict(st.verdict); }
}
async function enginePick() {
  const st = S.studio, node = cur(), btn = $("#pickBtn"); btn.disabled = true; btn.textContent = "choosing…";
  try {
    const r = await api("/api/studio/pick", { state: node.state, rule: st.rule, exclude: st.trail.map((n) => n.state.smiles) });
    if (S.studio !== st || cur() !== node) return;
    if (!r.ok) { toast(r.reason || "No move available"); btn.disabled = false; btn.textContent = "Let the engine choose"; return; }
    st.trail.push({ state: r.state, record: r.record, cat: r.catalogue }); st.site = null; st.verdict = null;
    renderStudio();
  } catch (e) { toast(e.message); btn.disabled = false; btn.textContent = "Let the engine choose"; }
}
