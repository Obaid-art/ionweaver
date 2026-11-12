"""The IonWeaver web API, and the page that uses it.

    GET  /api/bootstrap         what the page needs once: families, operators, checks, numbers
    GET  /api/family/{id}       the editable ions of one published family, drawn
    POST /api/identify          here is my ion: which family, which core?
    POST /api/start             a starting ion with its protected core and what each operator can offer
    POST /api/generate          N distinct certified products of one ion            (job)
    POST /api/batch             a quota per family over several families           (job)
    POST /api/chain             sequential editing, this deep, under this rule     (job)
    POST /api/sheet             my own family sheet: curate, discover cores, build  (job)
    GET  /api/job/{id}          progress and result
    POST /api/export            the products as CSV, SDF or SMILES
    POST /api/depict            one structure, drawn with the role colours

Every structure the engine hands back was certified before it was offered;
nothing here repairs, retries or filters.
"""
from __future__ import annotations

import base64
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from rdkit import Chem

from ionweaver_web import WEB_DIR, __version__, engine_fingerprint
from ionweaver_web.descriptors import BODY_RULES, DESCRIPTORS, ELEMENT_ORDER
from ionweaver_web.exports import products_csv, products_sdf, products_smi
from ionweaver_web.generate import (OPERATOR_INFO, OPERATORS, REQUIREMENTS, RULES, batch, chain, component_svg,
                                    generate, published_workspace, start_payload)
from ionweaver_web.identify import identify, make_record
from ionweaver_web.jobs import store
from ionweaver_web.library import SEED, load_published
from ionweaver_web.render import ROLES, molecule_svg
from ionweaver_web.sheet import NOVEL_EXAMPLE_CSV, TEMPLATE_CSV, get_sheet, load_sheet
from ionweaver_web.studio import catalogue, engine_pick, focus_svg, state_from_record, try_move

app = FastAPI(title="IonWeaver", version=__version__,
              description="Certified graph editing of individual ions: new ions that keep their family.")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

P = load_published()
WS = published_workspace(P)
if not P.checks["all"]:
    raise SystemExit(f"the rebuilt library does not reproduce the archive: {P.checks}")

_FAMILY_SVG: Dict[int, str] = {}


# -- request bodies -------------------------------------------------------------------

class Start(BaseModel):
    """Any of the ways of naming a starting ion."""
    kind: str = "published"                     # published | own | sheet | state
    inventory_id: Optional[int] = None          # published or sheet ion
    sheet_id: Optional[str] = None
    smiles: Optional[str] = Field(default=None, max_length=4000)
    family_id: Optional[str] = None
    family_label: Optional[str] = None
    grammar_id: Optional[str] = None
    pattern: Optional[str] = None
    core_atoms: Optional[List[int]] = None
    workspace: Optional[str] = None             # "published" or "sheet:<id>" (kind = state)


class StudioRequest(BaseModel):
    state: Dict[str, Any]
    move: Optional[Dict[str, Any]] = None
    operators: List[str] = Field(default_factory=lambda: list(OPERATORS))
    exclude: List[str] = Field(default_factory=list, max_length=5000)
    rule: str = "random"
    focus: List[int] = Field(default_factory=list)


class PieceRequest(BaseModel):
    keys: List[str] = Field(default_factory=list, max_length=600)


Constraints = Optional[Dict[str, List[Union[float, str, None]]]]   # bounds [lo, hi], or allowed element symbols


class GenerateRequest(BaseModel):
    start: Start
    want: int = Field(default=50, ge=1, le=2000)
    operators: List[str] = Field(default_factory=lambda: list(OPERATORS))
    per_operator: Optional[int] = Field(default=None, ge=1, le=500)
    constraints: Constraints = None            # {"mw": [null, 250], "hbd": [0, 1]}
    cursor: Optional[Dict[str, int]] = None    # continue an earlier request
    exclude: List[str] = Field(default_factory=list, max_length=20000)


class BatchRequest(BaseModel):
    family_ids: List[str] = Field(default_factory=list)   # empty means every published family
    sheet_id: Optional[str] = None
    per_family: int = Field(default=25, ge=1, le=1000)
    ions_per_family: int = Field(default=5, ge=1, le=100)
    operators: List[str] = Field(default_factory=lambda: list(OPERATORS))
    constraints: Constraints = None
    svg: bool = False


class ChainRequest(BaseModel):
    start: Start
    rule: str = "local"
    depth: int = Field(default=6, ge=1, le=12)
    allowance: int = Field(default=50, ge=4, le=400)
    operators: List[str] = Field(default_factory=lambda: list(OPERATORS))


class SheetRequest(BaseModel):
    filename: str = Field(default="", max_length=200)
    content_base64: Optional[str] = Field(default=None, max_length=9_000_000)   # ~6 MB of file
    text: Optional[str] = Field(default=None, max_length=6_000_000)
    example: Union[bool, str] = False          # True / "template", or "novel"


class ExportRequest(BaseModel):
    products: List[Dict[str, Any]] = Field(default_factory=list)
    fmt: str = "csv"
    name: str = "ionweaver"


class DepictRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=4000)
    core_atoms: List[int] = Field(default_factory=list)
    added_atoms: List[int] = Field(default_factory=list)
    removed_atoms: List[int] = Field(default_factory=list)
    q_labels: bool = False
    number_body: bool = False
    width: int = Field(default=420, ge=120, le=1200)
    height: int = Field(default=300, ge=90, le=900)


# -- resolving a start ---------------------------------------------------------------------

def _workspace_named(name: Optional[str]):
    if name and name.startswith("sheet:"):
        sheet = get_sheet(name.split(":", 1)[1])
        if sheet is None:
            raise HTTPException(404, "That uploaded sheet is no longer in memory. Upload it again.")
        return sheet.workspace
    return WS


def _resolve(start: Start):
    """A starting-ion record and the workspace its edits are decided in."""
    if start.kind == "state":
        # an ion the studio built: its core atoms are known from the certification labels
        if not (start.smiles and start.core_atoms and start.family_id and start.pattern):
            raise HTTPException(422, "A built ion needs its SMILES, core atoms, family and core pattern.")
        ws = _workspace_named(start.workspace)
        try:
            rec = make_record(start.smiles, start.family_id, start.family_label or start.family_id,
                              start.grammar_id or "G0", start.pattern, start.core_atoms,
                              inventory_id=start.inventory_id if start.inventory_id is not None else -1)
        except (ValueError, TypeError) as exc:
            raise HTTPException(422, str(exc)) from exc
        return rec, ws
    if start.kind == "sheet":
        sheet = get_sheet(start.sheet_id or "")
        if sheet is None:
            raise HTTPException(404, "That uploaded sheet is no longer in memory. Upload it again.")
        try:
            return sheet.record(int(start.inventory_id)), sheet.workspace
        except (KeyError, TypeError) as exc:
            raise HTTPException(422, str(exc)) from exc
    if start.kind == "own":
        if not start.smiles:
            raise HTTPException(422, "Paste a SMILES string.")
        verdict = identify(start.smiles, P)
        if not verdict.get("ok"):
            raise HTTPException(422, verdict.get("reason", "This ion cannot be edited."))
        choice = verdict["best"]
        if start.pattern:
            found = next((m for m in [verdict["best"]] + verdict["alternatives"]
                          if m["pattern"] == start.pattern), None)
            if found is None:
                raise HTTPException(422, "This ion does not carry the chosen core with all of its charge on it.")
            choice = found
        if verdict.get("published_id") is not None and choice["pattern"] == P.record(verdict["published_id"])["core_pattern"]:
            return P.record(verdict["published_id"]), WS
        rec = make_record(verdict["canonical_smiles"], choice["family_id"], choice["family_label"],
                          choice["grammar_id"], choice["pattern"], choice["core_atoms"], inventory_id=-1)
        return rec, WS
    try:
        return P.record(int(start.inventory_id)), WS
    except (KeyError, TypeError) as exc:
        raise HTTPException(404, str(exc)) from exc


# -- what the page needs once ----------------------------------------------------------------

_REP_SVG: Dict[int, str] = {}


def _families_with_pictures() -> List[Dict[str, Any]]:
    out = []
    for fam in P.families:
        row = dict(fam)
        rid = fam.get("representative_id")
        if rid is not None:
            if rid not in _REP_SVG:
                r = P.record(rid)
                _REP_SVG[rid] = molecule_svg(Chem.MolFromSmiles(r["canonical_smiles"]),
                                             core=[int(i) for i in str(r["core_atoms"]).split()],
                                             max_width=210, max_height=140)
            row["representative_svg"] = _REP_SVG[rid]
            row["representative_smiles"] = P.record(rid)["canonical_smiles"]
        out.append(row)
    return out


@app.get("/api/bootstrap")
def api_bootstrap() -> Dict[str, Any]:
    return {
        "version": __version__, "engine": engine_fingerprint(), "library_checks": P.checks,
        "load_seconds": P.load_seconds, "counts": P.counts(), "headline": P.headline,
        "sensitivity": P.sensitivity, "string_controls": P.string_controls,
        "families": _families_with_pictures(),
        "operators": [{"key": k, **v} for k, v in OPERATOR_INFO.items()],
        "requirements": REQUIREMENTS, "rules": [{"key": k, **v} for k, v in RULES.items()],
        "roles": ROLES, "descriptors": DESCRIPTORS, "body_rules": BODY_RULES,
        "elements": sorted({Chem.GetPeriodicTable().GetElementSymbol(int(z)) for z in P.elements}
                           | {a.GetSymbol() for r in P.starting for a in Chem.MolFromSmiles(r["canonical_smiles"]).GetAtoms()},
                          key=lambda e: (ELEMENT_ORDER.index(e) if e in ELEMENT_ORDER else 99, e)),
        "seed": SEED, "curation": P.curation_summary,
    }


@app.get("/api/family/{family_id}")
def api_family(family_id: str) -> Dict[str, Any]:
    fam = P.family(family_id)
    if fam is None:
        raise HTTPException(404, "No such family.")
    ions = []
    for r in P.ions_of(family_id):
        iid = r["inventory_id"]
        if iid not in _FAMILY_SVG:
            core = [int(i) for i in str(r["core_atoms"]).split()]
            _FAMILY_SVG[iid] = molecule_svg(Chem.MolFromSmiles(r["canonical_smiles"]), core=core,
                                            max_width=230, max_height=160)
        ions.append({"inventory_id": iid, "smiles": r["canonical_smiles"], "grammar_id": r["grammar_id"],
                     "heavy_atoms": r["heavy_atoms"], "body_heavy_atoms": r["body_heavy_atoms"],
                     "body_components": r["body_components"], "svg": _FAMILY_SVG[iid]})
    return {**fam, "ions": ions}


@app.get("/api/library")
def api_library() -> Dict[str, Any]:
    """The component library: every side group the curated ions display, drawn."""
    labels = {f["family_id"]: f["family_label"] for f in P.families}
    first: Dict[str, set] = {}
    for inst in P.instances:
        if inst["admitted"]:
            first.setdefault(inst["component_key"], set()).add(labels.get(inst["family_id"], inst["family_id"]))
    rows = []
    for e in P.library:
        rows.append({
            "library_id": e["library_id"], "key": e["component_key"], "heavy_atoms": e["heavy_atoms"],
            "attachment_order": e["attachment_bond_orders"], "rings": e["rings_in_component"],
            "occurrences": e["occurrences"], "ions": e["ions"], "families": e["families"],
            "family_labels": sorted(first.get(e["component_key"], ())),
            "graft": bool(e["graft_eligible"]), "swap": bool(e["swap_eligible"]),
            "prune": bool(e["prune_eligible"]), "crossover": bool(e["crossover_eligible"]),
            "reason": "" if e["graft_eligible"] else e.get("graft_reason", ""),
            "terminal_fragments": e["terminal_fragments_at_two_or_more"],
            "svg": component_svg(e["component_key"]),
        })
    return {"components": rows, "fragments": len(P.fragments), "summary": P.library_summary}


def _prewarm() -> None:
    """Draw every editable ion and library component once, so galleries open at once."""
    import threading

    def run() -> None:
        for r in P.starting:
            iid = r["inventory_id"]
            if iid not in _FAMILY_SVG:
                core = [int(i) for i in str(r["core_atoms"]).split()]
                _FAMILY_SVG[iid] = molecule_svg(Chem.MolFromSmiles(r["canonical_smiles"]), core=core,
                                                max_width=230, max_height=160)
        for e in P.library:
            component_svg(e["component_key"])

    threading.Thread(target=run, daemon=True, name="prewarm").start()


_prewarm()


# -- naming the ion --------------------------------------------------------------------------

@app.post("/api/identify")
def api_identify(start: Start) -> Dict[str, Any]:
    if not start.smiles:
        raise HTTPException(422, "Paste a SMILES string.")
    verdict = identify(start.smiles, P)
    if verdict.get("ok"):
        mol = Chem.MolFromSmiles(verdict["canonical_smiles"])
        verdict["svg"] = molecule_svg(mol, core=verdict["best"]["core_atoms"], q_labels=True, max_width=440, max_height=320)
        for alt in verdict["alternatives"]:
            alt["svg"] = molecule_svg(Chem.MolFromSmiles(verdict["canonical_smiles"]), core=alt["core_atoms"],
                                      max_width=220, max_height=160)
    return verdict


@app.post("/api/start")
def api_start(start: Start) -> Dict[str, Any]:
    rec, ws = _resolve(start)
    return start_payload(rec, ws)


# -- the three ways to edit -----------------------------------------------------------------------

@app.post("/api/generate")
def api_generate(req: GenerateRequest) -> Dict[str, Any]:
    rec, ws = _resolve(req.start)
    ops = [o for o in req.operators if o in OPERATORS] or list(OPERATORS)

    def work(progress):
        progress(0.1, "deciding every definition the operators offer")
        out = generate(rec, ws, want=req.want, operators=ops, per_operator=req.per_operator,
                       constraints=req.constraints, cursor=req.cursor, exclude=req.exclude)
        out["start"] = start_payload(rec, ws, width=420, height=300)
        return out

    return store.submit("generate", f"{req.want} ions", work).status()


@app.post("/api/batch")
def api_batch(req: BatchRequest) -> Dict[str, Any]:
    ops = [o for o in req.operators if o in OPERATORS] or list(OPERATORS)
    if req.sheet_id:
        sheet = get_sheet(req.sheet_id)
        if sheet is None:
            raise HTTPException(404, "That uploaded sheet is no longer in memory. Upload it again.")
        fids = req.family_ids or sorted({r["family_id"] for r in sheet.starting})
        ws, records_of = sheet.workspace, sheet.ions_of
    else:
        fids = req.family_ids or [f["family_id"] for f in P.families if f["editable"]]
        ws, records_of = WS, None
    if not fids:
        raise HTTPException(422, "No family has an editable ion.")

    def work(progress):
        return batch(P, ws, family_ids=fids, per_family=req.per_family, operators=ops,
                     ions_per_family=req.ions_per_family, records_of=records_of,
                     constraints=req.constraints, progress=progress, svg=req.svg)

    return store.submit("batch", f"{len(fids)} families × {req.per_family}", work).status()


@app.post("/api/chain")
def api_chain(req: ChainRequest) -> Dict[str, Any]:
    rec, ws = _resolve(req.start)
    if req.rule not in RULES:
        raise HTTPException(422, "Unknown selection rule.")
    ops = [o for o in req.operators if o in OPERATORS] or list(OPERATORS)

    def work(progress):
        out = chain(rec, ws, rule=req.rule, depth=req.depth, allowance=req.allowance,
                    operators=ops, progress=progress)
        out["start"] = start_payload(rec, ws, width=300, height=220)
        return out

    return store.submit("chain", f"{RULES[req.rule]['label']}, {req.depth} rounds", work).status()


# -- the studio: one chosen edit at a time ---------------------------------------------------------

@app.post("/api/studio/open")
def api_studio_open(start: Start) -> Dict[str, Any]:
    rec, ws = _resolve(start)
    return catalogue(state_from_record(rec, ws), ws)


def _studio(fn):
    """Run a studio call; a malformed state or move is a 422, never a 500."""
    try:
        return fn()
    except HTTPException:
        raise
    except (KeyError, ValueError, TypeError, IndexError) as exc:
        raise HTTPException(422, f"That request is not well formed: {exc}") from exc


@app.post("/api/studio/catalogue")
def api_studio_catalogue(req: StudioRequest) -> Dict[str, Any]:
    ws = _workspace_named(req.state.get("workspace"))
    return _studio(lambda: catalogue(req.state, ws))


@app.post("/api/studio/focus")
def api_studio_focus(req: StudioRequest) -> Dict[str, str]:
    ws = _workspace_named(req.state.get("workspace"))
    return _studio(lambda: {"svg": focus_svg(req.state, ws, req.focus)})


@app.post("/api/studio/try")
def api_studio_try(req: StudioRequest) -> Dict[str, Any]:
    if not req.move:
        raise HTTPException(422, "Name a move.")
    ws = _workspace_named(req.state.get("workspace"))

    def go():
        out = try_move(req.state, ws, req.move)
        if out.get("ok"):
            out["catalogue"] = catalogue(out["state"], ws)
        return out
    return _studio(go)


@app.post("/api/studio/pick")
def api_studio_pick(req: StudioRequest) -> Dict[str, Any]:
    ws = _workspace_named(req.state.get("workspace"))
    ops = [o for o in req.operators if o in OPERATORS] or list(OPERATORS)

    def go():
        out = engine_pick(req.state, ws, operators=ops, exclude=req.exclude, rule=req.rule)
        if out.get("ok"):
            out["catalogue"] = catalogue(out["state"], ws)
        return out
    return _studio(go)


@app.post("/api/piece")
def api_piece(req: PieceRequest) -> Dict[str, str]:
    """Library components or tails, drawn on their own."""
    return {k: component_svg(k) for k in req.keys[:600] if k}


# -- the user's own family ----------------------------------------------------------------------------

@app.post("/api/sheet")
def api_sheet(req: SheetRequest) -> Dict[str, Any]:
    if req.example == "novel":
        data, filename = NOVEL_EXAMPLE_CSV.encode("utf-8"), "families_nobody_published.csv"
    elif req.example:
        data, filename = TEMPLATE_CSV.encode("utf-8"), "example_family_sheet.csv"
    elif req.content_base64:
        try:
            data = base64.b64decode(req.content_base64)
        except Exception as exc:
            raise HTTPException(400, "The file could not be read.") from exc
        filename = req.filename or "uploaded.csv"
    elif req.text:
        data, filename = req.text.encode("utf-8"), req.filename or "pasted.csv"
    else:
        raise HTTPException(400, "Upload a file, paste rows, or try the example.")

    def work(progress):
        progress(0.1, "reading and curating every row")
        sheet = load_sheet(P, data, filename)
        progress(0.8, "drawing the cores")
        payload = sheet.as_dict()
        for fam in payload["families"]:
            for core in fam["cores"]:
                rep = sheet.record(core["copied_from"]) if core["copied_from"] in sheet.by_id else None
                if rep is not None:
                    mol = Chem.MolFromSmiles(rep["canonical_smiles"])
                    core["representative_smiles"] = rep["canonical_smiles"]
                    core["svg"] = molecule_svg(mol, core=[int(i) for i in str(rep["core_atoms"]).split()],
                                               q_labels=True, max_width=320, max_height=220)
        for ion in payload["ions"]:
            ion["svg"] = molecule_svg(Chem.MolFromSmiles(ion["canonical_smiles"]), core=ion["core_atoms"],
                                      max_width=230, max_height=160)
        return payload

    return store.submit("sheet", filename, work).status()


@app.get("/api/example_novel.csv")
def api_example_novel() -> Response:
    return Response(NOVEL_EXAMPLE_CSV, media_type="text/csv",
                    headers={"Content-Disposition": 'attachment; filename="families_nobody_published.csv"'})


@app.get("/api/template.csv")
def api_template() -> Response:
    return Response(TEMPLATE_CSV, media_type="text/csv",
                    headers={"Content-Disposition": 'attachment; filename="ionweaver_family_template.csv"'})


# -- jobs, files, pictures -----------------------------------------------------------------------------

@app.get("/api/job/{job_id}")
def api_job(job_id: str) -> Dict[str, Any]:
    job = store.get(job_id)
    if job is None:
        raise HTTPException(404, "Unknown job.")
    payload = job.status()
    if job.state == "done":
        payload["result"] = job.result
    return payload


@app.post("/api/export")
def api_export(req: ExportRequest) -> Response:
    name = "".join(c for c in req.name if c.isalnum() or c in "-_") or "ionweaver"
    if req.fmt == "sdf":
        return Response(products_sdf(req.products), media_type="chemical/x-mdl-sdfile",
                        headers={"Content-Disposition": f'attachment; filename="{name}.sdf"'})
    if req.fmt == "smi":
        return Response(products_smi(req.products), media_type="text/plain",
                        headers={"Content-Disposition": f'attachment; filename="{name}.smi"'})
    return Response(products_csv(req.products), media_type="text/csv",
                    headers={"Content-Disposition": f'attachment; filename="{name}.csv"'})


@app.post("/api/depict")
def api_depict(req: DepictRequest) -> Dict[str, str]:
    mol = Chem.MolFromSmiles(req.smiles)
    if mol is None:
        raise HTTPException(422, "That SMILES could not be read.")
    notes = ({i: str(i) for i in range(mol.GetNumAtoms()) if i not in set(req.core_atoms)}
             if req.number_body else None)
    return {"svg": molecule_svg(mol, core=req.core_atoms, added=req.added_atoms, removed=req.removed_atoms,
                                q_labels=req.q_labels, notes=notes, max_width=req.width, max_height=req.height)}


@app.get("/api/health")
def api_health() -> Dict[str, Any]:
    return {"ok": True, "version": __version__, "engine_unmodified": engine_fingerprint()["unmodified"]}


# -- the page -------------------------------------------------------------------------------------------

@app.get("/")
def index() -> FileResponse:
    return FileResponse(WEB_DIR / "index.html")


app.mount("/", StaticFiles(directory=str(WEB_DIR)), name="web")
