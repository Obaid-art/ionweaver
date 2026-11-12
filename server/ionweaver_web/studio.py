"""The studio: build an ion one chosen edit at a time.

A user points at a site on the current ion and asks for a particular move: a
library component at this atom, this element appended here, this side group
exchanged for that one.  The answer comes from the engine and nowhere else:

* the move is looked up among the definitions the engine's enumerators offer
  on this ion and the conditions admit (``operators.definitions_for``);
* if it is there, it is realised and certified (``build.realise``,
  ``build.certify``) and the product becomes the next ion, its core carried by
  the labels certification read;
* if it is not there, the same enumerator is run for that site and that
  structure alone and ``apply_conditions`` names the condition that refused it
  and why; a site no enumerator offers at all is explained from the rules that
  keep it out (a core atom, a ring bond, a one-atom component, an order that
  does not match).

Every product of the studio is therefore a product the one-edit generator
could also have offered; the studio only lets a person choose which.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem

from iw import build as B
from iw import operators as O
from iw import single as S
from iw.common import STUDY_SCORING, judge, parse
from iw.fragments import MIN_FRAGMENT_HEAVY_ATOMS, MIN_GROUP_HEAVY_ATOMS, terminal_pieces, _copy_out

from .generate import OPERATORS, Workspace, _record, capacity, component_svg
from .library import SEED
from .render import molecule_svg

__all__ = ["state_from_record", "catalogue", "try_move", "engine_pick", "move_key"]

PT = Chem.GetPeriodicTable()


# -- the state a client carries between calls ---------------------------------------------

def state_from_record(record: Dict[str, Any], ws: Workspace) -> Dict[str, Any]:
    return {
        "smiles": record["canonical_smiles"],
        "core_atoms": [int(i) for i in str(record["core_atoms"]).split()],
        "family_id": record["family_id"], "family_label": record.get("family_label", ""),
        "grammar_id": record["grammar_id"], "core_pattern": record["core_pattern"],
        "inventory_id": int(record["inventory_id"]), "workspace": ws.name,
        "origin_smiles": record["canonical_smiles"], "step": 0,
    }


def _context(state: Dict[str, Any], ws: Workspace) -> O.Context:
    for key in ("smiles", "core_atoms", "family_id", "grammar_id", "core_pattern"):
        if key not in state:
            raise ValueError(f"the ion state has no {key}")
    try:
        mol = parse(str(state["smiles"]))
    except Exception as exc:  # noqa: BLE001
        raise ValueError("the ion's SMILES could not be read") from exc
    core = sorted({int(i) for i in state["core_atoms"]})
    if not core or core[0] < 0 or core[-1] >= mol.GetNumAtoms():
        raise ValueError("the ion's core atoms do not index its atoms")
    if not str(state["family_id"]).startswith(("cation", "anion")):
        raise ValueError("the family must be a cation or an anion family")
    return O.context_from(mol, core, int(state.get("inventory_id", -1)),
                          str(state["family_id"]), Chem.MolToSmiles(mol), ws.elements)


def _definitions(ion: O.Context, ws: Workspace) -> List[O.Definition]:
    kept, _tally, _refused = O.definitions_for(ion, ws.library, ws.fragments, ws.observed, ws.observed_bonds)
    return kept


def _renumber(product_mol: Chem.Mol) -> Tuple[str, List[int]]:
    """The product as canonical SMILES, with its core atoms in that SMILES' atom order."""
    smiles = Chem.MolToSmiles(product_mol)
    order = list(product_mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    core_now = [int(a.GetIdx()) for a in product_mol.GetAtoms()
                if a.HasProp(B.CORE) and int(a.GetIntProp(B.CORE)) == 1]
    position = {old: new for new, old in enumerate(order)}
    return smiles, sorted(position[i] for i in core_now)


# -- naming a move -------------------------------------------------------------------------

def move_key(d: O.Definition) -> Tuple[Any, ...]:
    """A definition in the words a user chooses it by."""
    if d.operator == "Insert" and d.mode == "appended":
        (atom, order), = d.attachments[0].anchors
        return ("Insert", "appended", int(atom), PT.GetElementSymbol(d.builds[0].atomic_number), order)
    if d.operator == "Insert" and d.mode == "bridged":
        anchors = tuple(sorted(int(a) for a, _o in d.attachments[0].anchors))
        return ("Insert", "bridged", anchors, PT.GetElementSymbol(d.builds[0].atomic_number))
    if d.operator == "Substitute":
        atom, z = d.element_changes[0]
        return ("Substitute", int(atom), PT.GetElementSymbol(int(z)))
    if d.operator == "Rebond":
        a, b, order = d.order_changes[0]
        return ("Rebond", tuple(sorted((int(a), int(b)))), order)
    if d.operator == "Delete":
        return ("Delete", int(next(iter(d.remove_atoms))))
    if d.operator == "Graft":
        (atom, _order), = d.attachments[0].anchors
        return ("Graft", int(atom), d.attachments[0].structure)
    if d.operator == "Prune":
        a, b = next(iter(d.remove_bonds))
        anchor = int(d.touches[0].index)
        root = int(b if a == anchor else a)
        return ("Prune", anchor, root)
    if d.operator == "Swap":
        (anchor, _order), = d.attachments[0].anchors
        a, b = next(iter(d.remove_bonds))
        root = int(b if a == anchor else a)
        return ("Swap", int(anchor), root, d.attachments[0].structure)
    if d.operator == "Cross-over":
        (anchor, _order), = d.attachments[0].anchors
        a, b = next(iter(d.remove_bonds))
        root = int(b if a == anchor else a)
        return ("Cross-over", int(anchor), root, d.attachments[0].structure)
    return (d.operator, d.mode, d.site)


def key_of_request(move: Dict[str, Any]) -> Tuple[Any, ...]:
    op = move.get("operator")
    if op == "Insert" and move.get("mode", "appended") == "appended":
        return ("Insert", "appended", int(move["atom"]), str(move["element"]), str(move.get("order", "SINGLE")))
    if op == "Insert":
        return ("Insert", "bridged", tuple(sorted(int(x) for x in move["bond"])), str(move["element"]))
    if op == "Substitute":
        return ("Substitute", int(move["atom"]), str(move["element"]))
    if op == "Rebond":
        return ("Rebond", tuple(sorted(int(x) for x in move["bond"])), str(move["order"]))
    if op == "Delete":
        return ("Delete", int(move["atom"]))
    if op == "Graft":
        return ("Graft", int(move["atom"]), str(move["structure"]))
    if op == "Prune":
        return ("Prune", int(move["anchor"]), int(move["root"]))
    if op == "Swap":
        return ("Swap", int(move["anchor"]), int(move["root"]), str(move["structure"]))
    if op == "Cross-over":
        return ("Cross-over", int(move["anchor"]), int(move["root"]), str(move["structure"]))
    raise ValueError("unknown move")


# -- what the current ion offers ----------------------------------------------------------------

def catalogue(state: Dict[str, Any], ws: Workspace, *, svg: bool = True) -> Dict[str, Any]:
    ion = _context(state, ws)
    kept = _definitions(ion, ws)
    mol = ion.mol
    comp_of = {}
    for ci, c in enumerate(ion.components):
        for a in c.atoms:
            comp_of[int(a)] = ci
    core = sorted(ion.core)
    atoms = []
    for a in mol.GetAtoms():
        i = int(a.GetIdx())
        atoms.append({"idx": i, "symbol": a.GetSymbol(), "core": i in ion.core,
                      "q": core.index(i) if i in ion.core else None, "fixed": i in ion.fixed,
                      "ring": a.IsInRing(), "aromatic": a.GetIsAromatic(), "charge": int(a.GetFormalCharge()),
                      "degree": a.GetDegree(), "hs": int(a.GetTotalNumHs()), "component": comp_of.get(i)})
    bonds = [{"a": int(b.GetBeginAtomIdx()), "b": int(b.GetEndAtomIdx()), "order": str(b.GetBondType()),
              "editable": ion.editable_bond(b), "ring": b.IsInRing()} for b in mol.GetBonds()]
    components = [{"idx": ci, "key": c.key, "atoms": sorted(c.atoms), "root": c.root, "anchor": c.anchor,
                   "order": c.order, "heavy_atoms": len(c.atoms)} for ci, c in enumerate(ion.components)]

    moves: Dict[str, Any] = {"append": {}, "substitute": {}, "delete": [], "graft": {}, "insert": {},
                             "rebond": {}, "prune": [], "swap": {}, "crossover": {}}
    for d in kept:
        k = move_key(d)
        if k[0] == "Insert" and k[1] == "appended":
            moves["append"].setdefault(str(k[2]), []).append([k[3], k[4]])
        elif k[0] == "Insert":
            moves["insert"].setdefault("%d-%d" % k[2], []).append(k[3])
        elif k[0] == "Substitute":
            moves["substitute"].setdefault(str(k[1]), []).append(k[2])
        elif k[0] == "Rebond":
            moves["rebond"].setdefault("%d-%d" % k[1], []).append(k[2])
        elif k[0] == "Delete":
            moves["delete"].append(k[1])
        elif k[0] == "Graft":
            moves["graft"].setdefault(str(k[1]), []).append(k[2])
        elif k[0] == "Prune":
            moves["prune"].append({"anchor": k[1], "root": k[2], "atoms": len(d.remove_atoms),
                                   "description": d.description})
        elif k[0] == "Swap":
            moves["swap"].setdefault(str(comp_of.get(k[2], k[1])), []).append(k[3])
        elif k[0] == "Cross-over":
            ci = comp_of.get(k[2], None)
            slot = moves["crossover"].setdefault(str(ci), {})
            cut = slot.setdefault("%d-%d" % (k[1], k[2]), {"anchor": k[1], "root": k[2],
                                                            "removed_atoms": len(d.remove_atoms), "fragments": []})
            cut["fragments"].append(k[3])
    for key in ("append", "substitute", "graft", "insert", "rebond", "swap"):
        for slot in moves[key].values():
            slot.sort(key=lambda x: (str(x[0]) if isinstance(x, list) else str(x)))
    moves["delete"].sort()

    out = {
        "state": state, "atoms": atoms, "bonds": bonds, "components": components, "moves": moves,
        "capacity": capacity(kept), "total": len(kept),
        "library": [{"key": e["component_key"], "heavy_atoms": e["heavy_atoms"], "order": e["attachment_bond_orders"],
                     "graft": bool(e["graft_eligible"]), "swap": bool(e["swap_eligible"])} for e in ws.library],
        "fragments": [{"key": f["fragment_key"], "heavy_atoms": f["heavy_atoms"], "order": f["cut_bond_order"]}
                      for f in ws.fragments],
        "heavy_atoms": int(mol.GetNumAtoms()),
        "elements": [PT.GetElementSymbol(int(z)) for z in ws.elements],
        "orders": ["SINGLE", "DOUBLE", "TRIPLE"],
    }
    if svg:
        notes = {i: str(i) for i in range(mol.GetNumAtoms()) if i not in ion.core}
        out["svg"] = molecule_svg(mol, core=core, q_labels=True, notes=notes, max_width=560, max_height=380)
    return out


def focus_svg(state: Dict[str, Any], ws: Workspace, focus: Iterable[int]) -> str:
    ion = _context(state, ws)
    core = sorted(ion.core)
    notes = {i: str(i) for i in range(ion.mol.GetNumAtoms()) if i not in ion.core}
    return molecule_svg(ion.mol, core=core, q_labels=True, notes=notes, focus=list(focus),
                        max_width=560, max_height=380)


# -- trying one move -------------------------------------------------------------------------------

def _explain_missing(ion: O.Context, ws: Workspace, move: Dict[str, Any], key: Tuple[Any, ...]) -> Dict[str, Any]:
    """Why the engine does not offer this move here."""
    op = key[0]
    mol = ion.mol

    def atom_words(i: int) -> str:
        a = mol.GetAtomWithIdx(int(i))
        return "atom %d (%s%s)" % (i, a.GetSymbol(), ", in the protected core" if int(i) in ion.core else "")

    # 1. structural reasons: the enumerator never offers it at this site
    if op in ("Substitute", "Delete"):
        i = key[1]
        a = mol.GetAtomWithIdx(int(i))
        if int(i) in ion.core:
            return {"condition": "protected core", "reason": "%s belongs to the protected core; no operator changes a core atom." % atom_words(i)}
        if int(i) in ion.fixed:
            return {"condition": "ring held with the core", "reason": "%s lies in a ring the core covers in part, so the ring is held whole." % atom_words(i)}
        if a.IsInRing() and op == "Substitute":
            return {"condition": "ring atom", "reason": "Substitute never changes a ring atom: an element change inside a ring can aromatise it with no bond touched."}
        if a.IsInRing():
            return {"condition": "ring atom", "reason": "Delete removes an atom of degree one or two outside a ring only; this atom is in a ring."}
        if a.GetIsAromatic() or int(a.GetFormalCharge()) != 0:
            return {"condition": "aromatic or charged", "reason": "Only a neutral, non-aromatic body atom may be substituted or deleted."}
        if op == "Delete" and a.GetDegree() > 2:
            return {"condition": "degree", "reason": "Delete handles an atom carrying one or two bonds; this one carries %d." % a.GetDegree()}
        if op == "Delete" and a.GetDegree() == 2:
            orders = {str(b.GetBondType()) for b in a.GetBonds()}
            if len(orders) > 1:
                return {"condition": "no shared order", "reason": "Its two bonds differ in order, so there is no shared order at which to rejoin its neighbours."}
    if op == "Rebond" or (op == "Insert" and key[1] == "bridged"):
        a, b = key[2] if op == "Insert" else key[1]
        bond = mol.GetBondBetweenAtoms(int(a), int(b))
        if bond is None:
            return {"condition": "no such bond", "reason": "Atoms %d and %d are not bonded." % (a, b)}
        if bond.IsInRing():
            return {"condition": "ring bond", "reason": "No operator alters a ring bond."}
        if int(a) in ion.core and int(b) in ion.core:
            return {"condition": "protected core", "reason": "The bond between atoms %d and %d is inside the protected core." % (a, b)}
    if op == "Graft":
        entry = next((e for e in ws.library if e["component_key"] == key[2]), None)
        if entry is None:
            return {"condition": "not in the library", "reason": "%s is not a component the curated ions display; Graft attaches library components only." % key[2]}
        if not entry["graft_eligible"]:
            return {"condition": "one atom", "reason": "A single-atom addition is an Insert, not a Graft."}
    if op == "Swap":
        comp = next((c for c in ion.components if c.anchor == key[1] and c.root == key[2]), None)
        entry = next((e for e in ws.library if e["component_key"] == key[3]), None)
        if comp is None:
            return {"condition": "no side group", "reason": "No side group runs from atom %d through atom %d." % (key[1], key[2])}
        if len(comp.atoms) < MIN_GROUP_HEAVY_ATOMS:
            return {"condition": "one atom", "reason": "The side group at atom %d is a single atom; exchanging it is a Substitute, not a Swap." % key[1]}
        if entry is None:
            return {"condition": "not in the library", "reason": "%s is not a library component." % key[3]}
        if entry["attachment_bond_orders"] != comp.order:
            return {"condition": "bond order", "reason": "The side group is attached by a %s bond and %s attaches by a %s bond; Swap keeps the bond order." % (comp.order.lower(), key[3], entry["attachment_bond_orders"].lower())}
        if entry["component_key"] == comp.key:
            return {"condition": "same component", "reason": "That is the side group already there."}
    if op == "Cross-over":
        anchor, root = key[1], key[2]
        comp = next((c for c in ion.components if int(anchor) in c.atoms and int(root) in c.atoms), None)
        if comp is None:
            return {"condition": "not a side-group bond", "reason": "Cross-over cuts a bond inside a side group; this bond is not inside one."}
        cuts = terminal_pieces(mol, comp.atoms, frozenset({comp.root}))
        cut = next((c for c in cuts if {c.anchor, c.root} == {int(anchor), int(root)}), None)
        if cut is None:
            return {"condition": "no terminal piece", "reason": "No terminal piece of at least %d atoms lies beyond this bond, away from the core." % MIN_FRAGMENT_HEAVY_ATOMS}
        frag = next((f for f in ws.fragments if f["fragment_key"] == key[3]), None)
        if frag is None:
            return {"condition": "not a library tail", "reason": "%s is not a terminal piece of any library component." % key[3]}
        if frag["cut_bond_order"] != str(cut.order):
            return {"condition": "bond order", "reason": "The tail would be joined by a %s bond and this cut is a %s bond." % (frag["cut_bond_order"].lower(), str(cut.order).lower())}
    if op == "Prune":
        anchor, root = key[1], key[2]
        comp = next((c for c in ion.components if int(anchor) in c.atoms and int(root) in c.atoms), None)
        if comp is None:
            return {"condition": "not a side-group bond", "reason": "Prune cuts a bond inside a side group; this bond is not inside one."}
        cuts = terminal_pieces(mol, comp.atoms, frozenset({comp.root}))
        if any({c.anchor, c.root} == {int(anchor), int(root)} for c in cuts):
            return {"condition": "orientation", "reason": "A terminal piece lies beyond this bond on the other side; the move is offered as anchor %d, root %d." % (cuts[0].anchor, cuts[0].root)}
        return {"condition": "no terminal piece", "reason": "No terminal piece of at least %d atoms lies beyond this bond, away from the core." % MIN_FRAGMENT_HEAVY_ATOMS}
    if op == "Insert" and key[1] == "appended":
        z = PT.GetAtomicNumber(key[3]) if key[3] in [PT.GetElementSymbol(e) for e in ws.elements] else None
        if z is None:
            return {"condition": "element", "reason": "%s is not an element the inventory lets an edit place (%s)." % (key[3], ", ".join(PT.GetElementSymbol(e) for e in ws.elements))}
    if op == "Substitute":
        if key[2] not in [PT.GetElementSymbol(e) for e in ws.elements]:
            return {"condition": "element", "reason": "%s is not an element the inventory lets an edit place (%s)." % (key[2], ", ".join(PT.GetElementSymbol(e) for e in ws.elements))}

    # 2. the conditions: run the enumerator for this site and structure alone
    streams = []
    if op == "Graft":
        streams = [O.graft(ion, [e for e in ws.library if e["component_key"] == key[2]])]
    elif op == "Swap":
        streams = [O.swap(ion, [e for e in ws.library if e["component_key"] == key[3]])]
    elif op == "Cross-over":
        streams = [O.crossover(ion, [f for f in ws.fragments if f["fragment_key"] == key[3]])]
    elif op == "Insert":
        streams = [O.insert_appended(ion)] if key[1] == "appended" else [O.insert_bridged(ion)]
    elif op == "Substitute":
        streams = [O.substitute(ion)]
    elif op == "Rebond":
        streams = [O.rebond(ion)]
    elif op == "Delete":
        streams = [O.delete(ion)]
    elif op == "Prune":
        streams = [O.prune(ion)]
    for stream in streams:
        for d in stream:
            if move_key(d) != key:
                continue
            decided = O.apply_conditions(ion, d, ws.observed, ws.observed_bonds)
            if decided.condition_failed:
                blame = next((r for r in decided.results if r.index == decided.condition_atom), None)
                env = ""
                if blame is not None:
                    z, charge, aromatic, valence, hydrogens, orders = blame.environment
                    env = "%s, charge %+d%s, valence %d, %d H, bonds %s" % (
                        PT.GetElementSymbol(z), charge, ", aromatic" if aromatic else "", valence, hydrogens,
                        " ".join(o.lower() for o in orders))
                where = ("the new atom" if decided.condition_atom < 0 else "atom %d" % decided.condition_atom)
                return {"condition": decided.condition_failed.replace("_", " "),
                        "reason": decided.condition_reason, "at": where, "environment": env,
                        "description": d.description}
            return {"condition": "offered", "reason": "The engine offers this move; it was not applied because the request did not match it exactly."}
    return {"condition": "not offered", "reason": "No enumerator offers this move at this site."}


def try_move(state: Dict[str, Any], ws: Workspace, move: Dict[str, Any], *, svg: bool = True) -> Dict[str, Any]:
    ion = _context(state, ws)
    key = key_of_request(move)
    kept = _definitions(ion, ws)
    found = next((d for d in kept if move_key(d) == key), None)
    if found is None:
        why = _explain_missing(ion, ws, move, key)
        return {"ok": False, "move": move, **why}
    return _apply(state, ws, ion, found, svg=svg)


def _apply(state: Dict[str, Any], ws: Workspace, ion: O.Context, definition: O.Definition, *,
           svg: bool, position: int = 0) -> Dict[str, Any]:
    labelled = B.label(ion.mol, sorted(ion.core))
    expected = B.core_induced_bonds(labelled, frozenset(ion.core))
    required = 1 if state["family_id"].startswith("cation") else -1
    product = B.realise(labelled, definition)
    if not product.built:
        return {"ok": False, "condition": product.failure, "reason": product.reason, "description": definition.description}
    cert = B.certify(product.mol, labelled, required, expected, len(ion.core))
    if not cert.certified:
        return {"ok": False, "condition": cert.failed, "reason": cert.reason, "description": definition.description,
                "product_smiles": product.smiles}
    query = Chem.MolFromSmarts(state["core_pattern"]) if state.get("core_pattern") else None
    pseudo = {"inventory_id": int(state.get("inventory_id", -1)), "family_id": state["family_id"],
              "family_label": state.get("family_label", ""), "grammar_id": state["grammar_id"],
              "core_pattern": state["core_pattern"]}
    rec = _record(ion, pseudo, definition, position, product, cert, ion.mol, query, required, ws,
                  svg=svg, width=300, height=220)
    smiles, core_atoms = _renumber(product.mol)
    new_state = {**state, "smiles": smiles, "core_atoms": core_atoms, "step": int(state.get("step", 0)) + 1}
    rec["similarity_to_origin"] = S.similarity(parse(state["origin_smiles"]), product.mol) if state.get("origin_smiles") else None
    rec["identical_to_held"] = product.smiles == ion.canonical_smiles
    return {"ok": True, "record": rec, "state": new_state}


def engine_pick(state: Dict[str, Any], ws: Workspace, *, operators: Sequence[str] = OPERATORS,
                exclude: Iterable[str] = (), rule: str = "random", svg: bool = True) -> Dict[str, Any]:
    """Let the engine choose: the first certified new product in the study's order
    (``random``), or the most / least similar to the ion the trail started from."""
    if rule not in ("random", "local", "exploration"):
        raise ValueError(f"unknown rule {rule!r}")
    ion = _context(state, ws)
    kept = [d for d in _definitions(ion, ws) if d.operator in set(operators)]
    seen = set(exclude) | {ion.canonical_smiles}
    ordered = S.order_for(SEED + int(state.get("step", 0)), int(state.get("inventory_id", -1)),
                          "studio|%s" % state["grammar_id"], kept, 60 if rule != "random" else None)
    best = None
    origin = parse(state["origin_smiles"]) if state.get("origin_smiles") else ion.mol
    for position, d in enumerate(ordered, start=1):
        res = _apply(state, ws, ion, d, svg=False, position=position)
        if not res["ok"] or res["record"]["product_smiles"] in seen:
            continue
        if rule == "random":
            best = (d, position)
            break
        sim = res["record"]["similarity_to_origin"]
        score = -sim if rule == "local" else sim
        if best is None or score < best[2]:
            best = (d, position, score)
    if best is None:
        return {"ok": False, "condition": "exhausted", "reason": "No certified product is available that is new to this trail."}
    return _apply(state, ws, ion, best[0], svg=svg, position=best[1])
