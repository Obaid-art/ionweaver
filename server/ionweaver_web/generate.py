"""Asking the engine for new ions.

Three questions, each answered with the engine's own functions and nothing else:

``generate``
    From one starting ion, N distinct certified products, one edit each.  The
    candidates of each operator are taken in the engine's keyed order
    (``single.order_for``), which is the order the submitted run examined them
    in, so for a reference ion under the run seed the first products are the
    published ones.  Operators are interleaved so a request for N ions gets a
    mix, and every delivered product is distinct: distinct from the parent, from
    every other product of the request, and from anything the caller excludes.

``batch``
    The same, for several starting ions of several families, with a quota per
    family.

``chain``
    Sequential editing: each selected product becomes the next starting
    structure, under one of the three selection rules of the study.  The loop
    is written here so the depth can be chosen; every step inside it is the
    engine's (``context_from``, ``definitions_for``, ``label``, ``realise``,
    ``certify``, ``order_for``, ``similarity``, ``pick_key``).

An edit is never repaired, retried or filtered after the fact.  A definition
the conditions refuse is never built; a product certification refuses is
recorded as a refusal and never delivered.
"""
from __future__ import annotations

import functools
import math
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

from iw import build as B
from iw import operators as O
from iw import sequential as Q
from iw import single as S
from iw.common import STUDY_SCORING, judge, parse

from .descriptors import body_profile, descriptors_of
from .library import SEED, Published
from .render import molecule_svg

__all__ = ["OPERATORS", "OPERATOR_INFO", "MODE_LABELS", "REQUIREMENTS", "RULES",
           "Workspace", "published_workspace", "definitions_for", "capacity",
           "start_payload", "generate", "batch", "chain"]

Progress = Optional[Callable[[float, str], None]]

OPERATORS: Tuple[str, ...] = O.ATOM_LEVEL + O.GROUP_LEVEL

OPERATOR_INFO: Dict[str, Dict[str, str]] = {
    "Insert": {"level": "atom", "text": "Puts one atom into an acyclic bond, or on the end of a chain."},
    "Substitute": {"level": "atom", "text": "Changes the element of one body atom, keeping its bonds."},
    "Rebond": {"level": "atom", "text": "Changes the order of one acyclic bond outside the core."},
    "Delete": {"level": "atom", "text": "Removes one body atom of degree two and rejoins its neighbours, or one of degree one."},
    "Graft": {"level": "group", "text": "Attaches a library component at a free valence."},
    "Prune": {"level": "group", "text": "Cuts a terminal piece of a side group away."},
    "Swap": {"level": "group", "text": "Exchanges one whole side group for a library component."},
    "Cross-over": {"level": "group", "text": "Keeps the part of a side group nearest the core and splices on another's tail."},
}

MODE_LABELS: Dict[str, str] = {
    "Insert:bridged": "insert between two bonded atoms",
    "Insert:appended": "add at the end of a chain",
    "Substitute:element": "change an element",
    "Rebond:order": "change a bond order",
    "Delete:degree_two": "remove and reconnect",
    "Delete:degree_one": "remove a terminal atom",
    "Graft:attach": "attach a component",
    "Prune:cut": "cut a terminal piece",
    "Swap:exchange": "exchange a side group",
    "Cross-over:splice": "splice on another tail",
}

REQUIREMENTS: List[Dict[str, str]] = [
    {"key": "connected_closed_shell", "short": "One molecule",
     "text": "One connected, closed-shell molecule that sanitises."},
    {"key": "required_charge", "short": "Charge",
     "text": "The formal charge its family requires (+1 or −1)."},
    {"key": "core_atoms", "short": "Core atoms",
     "text": "Every core atom present as it was: element, charge, aromaticity."},
    {"key": "core_bonds", "short": "Core bonds",
     "text": "The induced bond set on the core atoms unchanged, nothing missing and nothing added."},
    {"key": "charge_on_core", "short": "Charge on core",
     "text": "All formal charge inside the core."},
]

RULES: Dict[str, Dict[str, str]] = {
    "local": {"label": "Close analogues", "text": "Each round selects the certified product most similar to the starting ion."},
    "random": {"label": "Random", "text": "Each round selects a certified product with no distance preference."},
    "exploration": {"label": "Exploration", "text": "Each round selects the certified product least similar to the starting ion."},
}


@dataclass
class Workspace:
    """What definitions are decided against: a library, the observed sets, the elements."""
    name: str
    library: List[Dict[str, Any]]
    fragments: List[Dict[str, Any]]
    observed: FrozenSet
    observed_bonds: FrozenSet
    elements: Tuple[int, ...]
    known_smiles: Set[str] = field(default_factory=set)


def published_workspace(P: Published) -> Workspace:
    return Workspace("published", P.library, P.fragments, P.observed, P.observed_bonds,
                     tuple(P.elements), P.smiles_set)


def merged_workspace(P: Published, name: str, library: List[Dict[str, Any]],
                     fragments: List[Dict[str, Any]], observed: Iterable, observed_bonds: Iterable,
                     elements: Iterable[int], known: Iterable[str]) -> Workspace:
    """The published workspace plus what a user's own sheet displays."""
    keys = {e["component_key"] for e in P.library}
    lib = list(P.library) + [e for e in library if e["component_key"] not in keys]
    fkeys = {f["fragment_key"] for f in P.fragments}
    frags = list(P.fragments) + [f for f in fragments if f["fragment_key"] not in fkeys]
    return Workspace(name, lib, frags, P.observed | frozenset(observed),
                     P.observed_bonds | frozenset(observed_bonds),
                     tuple(sorted(set(P.elements) | set(elements))), set(P.smiles_set) | set(known))


# -- one starting ion -----------------------------------------------------------------

def context_of(record: Dict[str, Any], ws: Workspace) -> O.Context:
    return O.context(record, ws.elements)


_DEF_CACHE: Dict[Tuple[str, str, str], Tuple[List[O.Definition], Dict[Any, Any]]] = {}


def definitions_for(record: Dict[str, Any], ws: Workspace) -> Tuple[O.Context, List[O.Definition], Dict[Any, Any]]:
    """Every definition the eight operators offer on this ion, decided once."""
    ion = context_of(record, ws)
    key = (ws.name, record["canonical_smiles"], str(record["core_atoms"]))
    hit = _DEF_CACHE.get(key)
    if hit is None:
        kept, tally, _refused = O.definitions_for(ion, ws.library, ws.fragments,
                                                   ws.observed, ws.observed_bonds)
        hit = (kept, tally)
        if len(_DEF_CACHE) > 400:
            _DEF_CACHE.clear()
        _DEF_CACHE[key] = hit
    kept, tally = hit
    return ion, kept, tally


def capacity(kept: Sequence[O.Definition]) -> Dict[str, int]:
    out = {op: 0 for op in OPERATORS}
    for d in kept:
        out[d.operator] += 1
    return out


def mode_capacity(kept: Sequence[O.Definition]) -> Dict[str, List[Dict[str, Any]]]:
    counts: Dict[Tuple[str, str], int] = {}
    for d in kept:
        counts[(d.operator, d.mode)] = counts.get((d.operator, d.mode), 0) + 1
    out: Dict[str, List[Dict[str, Any]]] = {op: [] for op in OPERATORS}
    for (op, mode), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        out[op].append({"mode": mode, "label": MODE_LABELS.get(f"{op}:{mode}", mode), "count": n})
    return out


def side_groups(ion: O.Context, *, svg: bool = False) -> List[Dict[str, Any]]:
    out = []
    for c in ion.components:
        row = {"key": c.key, "atoms": sorted(c.atoms), "anchor_core_atom": c.anchor,
               "root": c.root, "order": c.order, "heavy_atoms": len(c.atoms)}
        if svg:
            row["svg"] = component_svg(c.key)
        out.append(row)
    return out


@functools.lru_cache(maxsize=4096)
def component_svg(key: str, width: int = 170, height: int = 120) -> str:
    """A library component or side group on its own, the attachment as a dummy."""
    mol = Chem.MolFromSmiles(key)
    if mol is None:
        return ""
    return molecule_svg(mol, max_width=width, max_height=height)


def start_payload(record: Dict[str, Any], ws: Workspace, *, svg: bool = True,
                  width: int = 440, height: int = 320) -> Dict[str, Any]:
    """The starting ion, its protected core, and what each operator can offer."""
    ion, kept, tally = definitions_for(record, ws)
    caps = capacity(kept)
    core = sorted(ion.core)
    payload = {
        "inventory_id": record["inventory_id"], "family_id": record["family_id"],
        "family_label": record.get("family_label", ""), "smiles": record["canonical_smiles"],
        "ion_type": "cation" if record["family_id"].startswith("cation") else "anion",
        "charge": 1 if record["family_id"].startswith("cation") else -1,
        "grammar_id": record["grammar_id"], "core_pattern": record["core_pattern"],
        "core_atoms": core, "ring_shared_atoms": sorted(ion.fixed - ion.core),
        "heavy_atoms": int(ion.mol.GetNumAtoms()), "body_heavy_atoms": int(ion.mol.GetNumAtoms() - len(ion.fixed)),
        "side_groups": side_groups(ion, svg=svg),
        "capacity": caps, "total_capacity": int(sum(caps.values())),
        "modes": mode_capacity(kept),
        "offered": {f"{op}:{mode}": row["offered"] for (op, mode), row in tally.items()},
        "descriptors": {**descriptors_of(ion.mol), **body_profile(ion.mol, core)},
        "workspace": ws.name,
    }
    if svg:
        payload["svg"] = molecule_svg(ion.mol, core=core, q_labels=True, max_width=width, max_height=height)
        payload["core_svg"] = core_only_svg(record["core_pattern"], ion.mol, core)
    return payload


@functools.lru_cache(maxsize=4096)
def _core_only(smiles: str, core_key: str) -> str:
    from rdkit.Chem import rdmolops
    mol = Chem.MolFromSmiles(smiles)
    core = [int(i) for i in core_key.split()]
    keep = set(core)
    rw = Chem.RWMol(mol)
    # mark the attachment points, then cut the body away
    for i in core:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            if int(nb.GetIdx()) not in keep:
                rw.GetAtomWithIdx(int(nb.GetIdx())).SetAtomicNum(0)
    for i in sorted((int(a.GetIdx()) for a in mol.GetAtoms() if int(a.GetIdx()) not in keep
                     and not any(int(n.GetIdx()) in keep for n in a.GetNeighbors())), reverse=True):
        rw.RemoveAtom(i)
    piece = rw.GetMol()
    try:
        piece.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(piece, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:  # noqa: BLE001 - fall back on the full parent
        return molecule_svg(mol, core=core, show_body=False, max_width=170, max_height=120)
    core_now = [int(a.GetIdx()) for a in piece.GetAtoms() if a.GetAtomicNum() != 0]
    return molecule_svg(piece, core=core_now, show_body=False, max_width=170, max_height=120)


def core_only_svg(pattern: str, mol: Chem.Mol, core: Sequence[int]) -> str:
    return _core_only(Chem.MolToSmiles(mol), " ".join(str(i) for i in core))


# -- roles of a product's atoms --------------------------------------------------------

def _roles(product_mol: Chem.Mol, ion: O.Context, definition: O.Definition) -> Dict[str, List[int]]:
    """Which product atoms are core, which the edit added, which it changed."""
    core, added, changed = [], [], []
    element_changed = {int(i) for i, _z in definition.element_changes}
    rebonded = {int(a) for a, _b, _o in definition.order_changes} | {int(b) for _a, b, _o in definition.order_changes}
    new_bonded = {int(a) for a, _b, _o in definition.add_bonds} | {int(b) for _a, b, _o in definition.add_bonds}
    for atom in product_mol.GetAtoms():
        idx = int(atom.GetIdx())
        if atom.HasProp(B.CORE) and int(atom.GetIntProp(B.CORE)) == 1:
            core.append(idx)
        elif atom.HasProp(B.ADDED):
            added.append(idx)
        elif atom.HasProp(B.PARENT):
            p = int(atom.GetIntProp(B.PARENT))
            if p in element_changed or p in rebonded or p in new_bonded:
                changed.append(idx)
        else:
            added.append(idx)
    return {"core": core, "added": added, "changed": changed}


def _canonical_positions(mol: Chem.Mol) -> Dict[int, int]:
    """Atom index in this Mol -> atom index in its canonical SMILES."""
    Chem.MolToSmiles(mol)
    order = list(mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    return {old: new for new, old in enumerate(order)}


def _canonical_core(product_mol: Chem.Mol) -> List[int]:
    """The product's core atoms, as indices into its canonical SMILES."""
    Chem.MolToSmiles(product_mol)
    order = list(product_mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    position = {old: new for new, old in enumerate(order)}
    return sorted(position[int(a.GetIdx())] for a in product_mol.GetAtoms()
                  if a.HasProp(B.CORE) and int(a.GetIntProp(B.CORE)) == 1)


def _cert_checks(cert: B.Certificate) -> Dict[str, bool]:
    """The five requirements as booleans, from what certification found."""
    reqs = [r["key"] for r in REQUIREMENTS]
    if cert.certified:
        return {k: True for k in reqs}
    failed_at = next((i for i, r in enumerate(B.REQUIREMENTS) if r == cert.failed), len(reqs))
    return {k: i < failed_at for i, k in enumerate(reqs)}


def _record(ion: O.Context, record: Dict[str, Any], definition: O.Definition, position: int,
            product: B.Product, cert: B.Certificate, parent_mol: Chem.Mol, query: Any,
            required: int, ws: Workspace, *, svg: bool, width: int, height: int) -> Dict[str, Any]:
    roles = _roles(product.mol, ion, definition)
    verdict = judge(product.smiles, ion.canonical_smiles, query, required, STUDY_SCORING)
    # the parent is reported as canonical SMILES, so parent atom indices are reported
    # in that SMILES' order (a chain's held ion is in build order, and differs)
    pos = _canonical_positions(ion.mol)
    out = {
        "inventory_id": record["inventory_id"], "family_id": record["family_id"],
        "family_label": record.get("family_label", ""), "grammar_id": record["grammar_id"],
        "core_pattern": record["core_pattern"], "parent_smiles": ion.canonical_smiles,
        "operator": definition.operator, "mode": definition.mode,
        "edit": MODE_LABELS.get(f"{definition.operator}:{definition.mode}", definition.mode),
        "site": definition.site, "description": definition.description,
        "signature": definition.signature_text, "order_index": position,
        "definition_id": "%d-%s-%d" % (ion.inventory_id, definition.operator, position),
        "heavy_atom_change": int(definition.heavy_atom_change),
        "product_smiles": product.smiles, "mapped_smiles": product.mapped_smiles,
        "similarity_to_parent": S.similarity(parent_mol, product.mol),
        "certified": bool(cert.certified), "failed": cert.failed, "reason": cert.reason,
        "checks": _cert_checks(cert),
        "judge_accepted": bool(verdict.accepted), "judge_reason": verdict.reason,
        "judge_agrees": bool(cert.certified) == bool(verdict.core_retained and verdict.charge_on_core
                                                     and verdict.required_charge and verdict.connected
                                                     and verdict.closed_shell),
        "new_to_inventory": product.smiles not in ws.known_smiles,
        "core_atoms": roles["core"], "added_atoms": roles["added"], "changed_atoms": roles["changed"],
        "parent_core_atoms": sorted(pos[int(i)] for i in ion.core),
        "product_core_atoms": _canonical_core(product.mol),
        "workspace": ws.name,
        "removed_parent_atoms": sorted(pos[int(i)] for i in definition.remove_atoms),
        "changed_parent_atoms": sorted({pos[int(i)] for i, _z in definition.element_changes}
                                       | {pos[int(x)] for a, b, _o in definition.order_changes for x in (a, b)}
                                       | {pos[int(x)] for a, b, _o in definition.add_bonds for x in (a, b)}
                                       | {pos[int(a)] for att in definition.attachments for a, _o in att.anchors}),
        "descriptors": {**descriptors_of(product.mol), **body_profile(product.mol, roles["core"])},
    }
    if svg:
        out["svg"] = molecule_svg(product.mol, core=roles["core"], added=roles["added"] + roles["changed"],
                                  max_width=width, max_height=height)
    return out


def parent_svg(record: Dict[str, Any], ws: Workspace, removed: Sequence[int] = (),
               changed: Sequence[int] = (), *, width: int = 300, height: int = 220) -> str:
    ion = context_of(record, ws)
    return molecule_svg(ion.mol, core=sorted(ion.core), removed=removed, added=changed,
                        max_width=width, max_height=height)


# -- N distinct products from one ion ----------------------------------------------------

Constraints = Optional[Dict[str, Sequence[Optional[float]]]]


def within(descriptors: Dict[str, Any], constraints: Constraints) -> bool:
    """Does a product's descriptor row sit inside every [low, high] asked for?"""
    if not constraints:
        return True
    for key, bounds in constraints.items():
        value = descriptors.get(key)
        if value is None:
            continue
        if key == "body_elements":
            allowed = {str(e) for e in (bounds or [])}
            if any(e not in allowed for e in value):
                return False
            continue
        lo, hi = (list(bounds) + [None, None])[:2]
        if lo is not None and value < float(lo) - 1e-9:
            return False
        if hi is not None and value > float(hi) + 1e-9:
            return False
    return True


@functools.lru_cache(maxsize=8192)
def _elements_in(structure: str) -> frozenset:
    """The element symbols a fragment or tail would bring in (the dummy excluded)."""
    mol = Chem.MolFromSmiles(structure, sanitize=False)
    if mol is None:
        return frozenset()
    return frozenset(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def introduces(definition: O.Definition) -> frozenset:
    """Every element an edit would place: a built atom, a new element, a joined structure."""
    out = {Chem.GetPeriodicTable().GetElementSymbol(int(b.atomic_number)) for b in definition.builds}
    out |= {Chem.GetPeriodicTable().GetElementSymbol(int(z)) for _i, z in definition.element_changes}
    for att in definition.attachments:
        out |= _elements_in(att.structure)
    return frozenset(out)


@functools.lru_cache(maxsize=8192)
def _structure_traits(structure: str) -> Tuple[bool, bool, bool]:
    """(has a double or triple bond, has a ring, has an aromatic ring) for a fragment or tail."""
    mol = Chem.MolFromSmiles(structure, sanitize=False)
    if mol is None:
        return (False, False, False)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    multiple = any(b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE) for b in mol.GetBonds())
    rings = mol.GetRingInfo().NumRings() > 0
    aromatic = any(a.GetIsAromatic() for a in mol.GetAtoms())
    return (multiple, rings, aromatic)


def _adds_multiple_bond(d: O.Definition) -> bool:
    if any(o in ("DOUBLE", "TRIPLE") for b in d.builds for o in b.bonds):
        return True
    if any(o in ("DOUBLE", "TRIPLE") for _a, _b, o in d.order_changes):
        return True
    if any(o in ("DOUBLE", "TRIPLE") for _a, _b, o in d.add_bonds):
        return True
    return any(_structure_traits(att.structure)[0] or any(o in ("DOUBLE", "TRIPLE") for _i, o in att.anchors)
               for att in d.attachments)


def prune_for_rules(kept: Sequence[O.Definition], constraints: Constraints, ion: O.Context) -> List[O.Definition]:
    """The definitions that can still obey the rules, decided before anything is built.

    A subset of what the engine offers, never an addition.  Only rules a
    definition decides on its own are applied here (a banned element, a new
    double or triple bond, a new ring, a size cap); every product is still
    checked in full after it is built.
    """
    rules = constraints or {}
    if not rules:
        return list(kept)
    hi = lambda k: (list(rules.get(k) or [None, None]) + [None, None])[1]
    allowed = {str(e) for e in rules["body_elements"]} if rules.get("body_elements") is not None else None
    no_multiple = hi("body_multiple_bonds") == 0
    no_rings = hi("body_rings") == 0
    no_aromatic = hi("body_aromatic_rings") == 0 or no_rings
    body_cap = hi("body_heavy_atoms")
    atoms_cap = hi("heavy_atoms")
    body_now = ion.mol.GetNumAtoms() - len(ion.fixed)
    atoms_now = ion.mol.GetNumAtoms()
    out = []
    for d in kept:
        if allowed is not None and not introduces(d) <= allowed:
            continue
        if no_multiple and _adds_multiple_bond(d):
            continue
        if no_rings and any(_structure_traits(att.structure)[1] for att in d.attachments):
            continue
        if no_aromatic and any(_structure_traits(att.structure)[2] for att in d.attachments):
            continue
        if body_cap is not None and body_now + int(d.heavy_atom_change) > float(body_cap):
            continue
        if atoms_cap is not None and atoms_now + int(d.heavy_atom_change) > float(atoms_cap):
            continue
        out.append(d)
    return out


def allowed_by_elements(kept: Sequence[O.Definition], allowed: Optional[Iterable[str]]) -> List[O.Definition]:
    """Kept for callers that only need the element rule."""
    if allowed is None:
        return list(kept)
    ok = {str(e) for e in allowed}
    return [d for d in kept if introduces(d) <= ok]


def generate(record: Dict[str, Any], ws: Workspace, *, want: int = 50,
             operators: Sequence[str] = OPERATORS, seed: int = SEED,
             per_operator: Optional[int] = None, exclude: Iterable[str] = (),
             constraints: Constraints = None, cursor: Optional[Dict[str, int]] = None,
             svg: bool = True, width: int = 300, height: int = 220,
             max_examinations: Optional[int] = None) -> Dict[str, Any]:
    """Up to ``want`` distinct certified products of one starting ion.

    Each operator's candidates are examined in the engine's keyed order, the
    order the submitted study examined them in; the operators are interleaved
    one candidate at a time.  A product is delivered only if certification
    passed, it is not the parent, not already delivered, not in ``exclude``,
    and inside ``constraints`` (descriptor bounds) when any are given.

    ``cursor`` continues an earlier request: each operator starts where the
    earlier one stopped, so "more ions" walks further down the same order.
    """
    t0 = time.perf_counter()
    ion, kept, _tally = definitions_for(record, ws)
    ops = [o for o in OPERATORS if o in set(operators)]
    usable = prune_for_rules(kept, constraints, ion)
    by_op = {op: [d for d in usable if d.operator == op] for op in ops}
    orders = {op: S.order_for(seed, ion.inventory_id, op, by_op[op], None) for op in ops}
    cursor = {op: int((cursor or {}).get(op, 0)) for op in ops}
    delivered_by_op = {op: 0 for op in ops}
    outside = 0
    set_aside = len(kept) - len(usable)

    labelled = B.label(ion.mol, sorted(ion.core))
    expected = B.core_induced_bonds(labelled, frozenset(ion.core))
    required = 1 if record["family_id"].startswith("cation") else -1
    query = Chem.MolFromSmarts(record["core_pattern"]) if record.get("core_pattern") else None
    parent_mol = parse(ion.canonical_smiles)

    seen: Set[str] = set(exclude) | {ion.canonical_smiles}
    products: List[Dict[str, Any]] = []
    refusals: List[Dict[str, Any]] = []
    examinations = 0
    repeats = 0
    identical = 0
    # the examination budget: generous under rules, since a bound on weight or donors
    # can only be decided after a product is built; a request that stops on budget
    # says so, and "more" continues from the same cursors
    limit = max_examinations if max_examinations is not None else max(400, (150 if constraints else 12) * int(want))
    want = max(0, int(want))

    live = [op for op in ops if orders[op]]
    while live and len(products) < want and examinations < limit:
        for op in list(live):
            if len(products) >= want or examinations >= limit:
                break
            if per_operator is not None and delivered_by_op[op] >= per_operator:
                live.remove(op)
                continue
            if cursor[op] >= len(orders[op]):
                live.remove(op)
                continue
            position = cursor[op] + 1
            definition = orders[op][cursor[op]]
            cursor[op] += 1
            examinations += 1
            product = B.realise(labelled, definition)
            if not product.built:
                refusals.append({"operator": op, "mode": definition.mode, "description": definition.description,
                                 "requirement": product.failure, "reason": product.reason})
                continue
            cert = B.certify(product.mol, labelled, required, expected, len(ion.core))
            if not cert.certified:
                refusals.append({"operator": op, "mode": definition.mode, "description": definition.description,
                                 "requirement": cert.failed, "reason": cert.reason,
                                 "product_smiles": product.smiles})
                continue
            if product.smiles == ion.canonical_smiles:
                identical += 1
                continue
            if product.smiles in seen:
                repeats += 1
                continue
            seen.add(product.smiles)
            rec = _record(ion, record, definition, position, product, cert, parent_mol, query,
                          required, ws, svg=svg, width=width, height=height)
            if not within(rec["descriptors"], constraints):
                outside += 1
                continue
            delivered_by_op[op] += 1
            products.append(rec)

    caps = capacity(kept)
    return {
        "start": {"inventory_id": record["inventory_id"], "smiles": ion.canonical_smiles,
                  "family_id": record["family_id"], "family_label": record.get("family_label", ""),
                  "grammar_id": record["grammar_id"], "core_pattern": record["core_pattern"],
                  "core_atoms": sorted(ion.core)},
        "products": products,
        "summary": {
            "requested": want, "delivered": len(products), "examinations": examinations,
            "refused_by_certification": len(refusals), "repeats": repeats, "identical_to_parent": identical,
            "outside_constraints": outside, "definitions_set_aside_by_rules": set_aside, "definitions_set_aside_by_element_rule": set_aside,
            "capacity": caps, "usable": {op: len(by_op[op]) for op in ops}, "examined_per_operator": dict(cursor),
            "cursor": dict(cursor),
            "delivered_per_operator": delivered_by_op,
            "exhausted": not live, "stopped_on_budget": bool(live) and len(products) < want and examinations >= limit,
            "budget": int(limit), "elapsed_s": round(time.perf_counter() - t0, 2),
            "new_to_inventory": sum(1 for p in products if p["new_to_inventory"]),
        },
        "refusals": refusals,
        "settings": {"seed": int(seed), "operators": ops, "per_operator": per_operator,
                     "constraints": dict(constraints or {}), "workspace": ws.name},
    }


# -- several ions of several families -----------------------------------------------------

def spread_by_size(ions: Sequence[Dict[str, Any]], n: int) -> List[Dict[str, Any]]:
    """``n`` ions of a family, spread evenly from its smallest body to its largest.

    Deterministic: the editable members are sorted by body size, then by
    canonical SMILES, and taken at evenly spaced positions, so a batch always
    uses the same parents and covers the family's range of sizes.
    """
    ordered = sorted(ions, key=lambda r: (r["body_heavy_atoms"], r["canonical_smiles"]))
    n = max(1, min(int(n), len(ordered)))
    if n == len(ordered):
        return ordered
    picks = sorted({int(round(i * (len(ordered) - 1) / max(1, n - 1))) for i in range(n)})
    return [ordered[i] for i in picks]


def batch(P: Published, ws: Workspace, *, family_ids: Sequence[str], per_family: int = 25,
          operators: Sequence[str] = OPERATORS, seed: int = SEED, ions_per_family: int = 5,
          records_of: Optional[Callable[[str], List[Dict[str, Any]]]] = None,
          constraints: Constraints = None,
          progress: Progress = None, svg: bool = False) -> Dict[str, Any]:
    """A quota of distinct products per family, spread over several of its ions.

    The parents of a family are ``ions_per_family`` of its editable members,
    spread evenly from the smallest body to the largest.  Products are gathered
    round-robin across those parents, so no single one dominates, and
    distinctness holds across the whole batch.
    """
    t0 = time.perf_counter()
    products: List[Dict[str, Any]] = []
    per_family_rows: List[Dict[str, Any]] = []
    seen: Set[str] = set()
    fam_ids = list(family_ids)
    for i, fid in enumerate(fam_ids):
        ions = list(records_of(fid)) if records_of else P.ions_of(fid)
        label = ions[0].get("family_label", fid) if ions else fid
        if progress:
            progress(i / max(1, len(fam_ids)), f"{label}: choosing starting ions")
        if not ions:
            per_family_rows.append({"family_id": fid, "family_label": label, "ions_used": 0,
                                    "delivered": 0, "requested": per_family, "note": "no editable ion"})
            continue
        chosen = spread_by_size(ions, ions_per_family)
        # parents whose own body already breaks a body rule: only an edit that removes
        # the offending atoms can comply, so such a family is likely to fall short
        body_rules = {k: v for k, v in (constraints or {}).items() if k.startswith("body_")}
        breaking = 0
        if body_rules:
            fine, poor = [], []
            for rec in chosen:
                prof = body_profile(parse(rec["canonical_smiles"]), [int(i) for i in str(rec["core_atoms"]).split()])
                (fine if within(prof, body_rules) else poor).append(rec)
            breaking = len(poor)
            # parents that already comply go first; the others are used only if the quota is still short
            chosen = fine + poor
        poor_ids = {id(r) for r in (poor if body_rules else [])}
        quota = int(per_family)
        share = int(math.ceil(quota / len(chosen)))
        got: List[Dict[str, Any]] = []
        cursors: Dict[int, Dict[str, int]] = {}
        exhausted: set = set()
        for j, rec in enumerate(chosen):
            if len(got) >= quota:
                break
            if progress:
                progress((i + j / len(chosen)) / max(1, len(fam_ids)),
                         f"{label}: ion {j + 1} of {len(chosen)}")
            # a parent whose body already breaks a rule rarely yields; give it a short look only
            res = generate(rec, ws, want=min(share, quota - len(got)), operators=operators, seed=seed,
                           exclude=seen, constraints=constraints, svg=svg,
                           max_examinations=600 if id(rec) in poor_ids else None)
            cursors[j] = res["summary"]["cursor"]
            if res["summary"]["exhausted"]:
                exhausted.add(j)
            for p in res["products"]:
                seen.add(p["product_smiles"])
                got.append(p)
        # further passes over the same parents, each continuing where it stopped,
        # until the quota is met or every parent is exhausted
        while len(got) < quota and len(exhausted) < len(chosen):
            before = len(got)
            for j, rec in enumerate(chosen):
                if len(got) >= quota or j in exhausted:
                    continue
                res = generate(rec, ws, want=quota - len(got), operators=operators, seed=seed,
                               exclude=seen, constraints=constraints, cursor=cursors.get(j), svg=svg,
                               max_examinations=600 if id(rec) in poor_ids else None)
                cursors[j] = res["summary"]["cursor"]
                if res["summary"]["exhausted"] or not res["products"]:
                    exhausted.add(j)
                for p in res["products"]:
                    seen.add(p["product_smiles"])
                    got.append(p)
            if len(got) == before:
                break
        products.extend(got)
        per_family_rows.append({
            "family_id": fid, "family_label": label, "ion_type": "cation" if fid.startswith("cation") else "anion",
            "ions_used": len(chosen), "ions_available": len(ions), "requested": quota, "delivered": len(got),
            "new_to_inventory": sum(1 for p in got if p["new_to_inventory"]),
            "operators": {op: sum(1 for p in got if p["operator"] == op) for op in OPERATORS},
            "parents_breaking_rules": breaking,
        })
    if progress:
        progress(1.0, "done")
    return {"products": products, "families": per_family_rows,
            "summary": {"families": len(fam_ids), "delivered": len(products),
                        "requested": int(per_family) * len(fam_ids),
                        "distinct": len({p["product_smiles"] for p in products}),
                        "new_to_inventory": sum(1 for p in products if p["new_to_inventory"]),
                        "elapsed_s": round(time.perf_counter() - t0, 2)},
            "settings": {"seed": int(seed), "per_family": int(per_family), "ions_per_family": int(ions_per_family),
                         "operators": list(operators), "constraints": dict(constraints or {}),
                         "workspace": ws.name}}


# -- sequential editing ------------------------------------------------------------------

def chain(record: Dict[str, Any], ws: Workspace, *, rule: str = "local", depth: int = 6,
          seed: int = SEED, allowance: int = 50, operators: Sequence[str] = OPERATORS,
          svg: bool = True, width: int = 260, height: int = 190, progress: Progress = None) -> Dict[str, Any]:
    """``depth`` rounds of sequential editing from one ion under one rule.

    This is the study's ``sequential.walk`` with the number of rounds and the
    operators as parameters.  Every chemistry call is the engine's.
    """
    if rule not in RULES:
        raise ValueError(f"unknown rule {rule!r}")
    t0 = time.perf_counter()
    preference = Q.PREFERENCE[rule]
    required = 1 if record["family_id"].startswith("cation") else -1
    held = parse(record["canonical_smiles"])
    start = parse(record["canonical_smiles"])
    held_core = sorted(int(i) for i in str(record["core_atoms"]).split())
    ops = set(operators)
    query = Chem.MolFromSmarts(record["core_pattern"]) if record.get("core_pattern") else None
    steps: List[Dict[str, Any]] = []
    seen: Set[str] = {Chem.MolToSmiles(held)}
    depth = max(1, min(int(depth), 12))

    for round_number in range(1, depth + 1):
        if progress:
            progress((round_number - 1) / depth, f"{RULES[rule]['label'].lower()}: round {round_number} of {depth}")
        held_smiles = Chem.MolToSmiles(held)
        ion = O.context_from(held, held_core, record["inventory_id"], record["family_id"], held_smiles, ws.elements)
        kept, _tally, _refused = O.definitions_for(ion, ws.library, ws.fragments, ws.observed, ws.observed_bonds)
        kept = [d for d in kept if d.operator in ops]
        labelled = B.label(ion.mol, sorted(ion.core))
        expected = B.core_induced_bonds(labelled, frozenset(ion.core))
        chosen = S.order_for(seed + round_number, record["inventory_id"],
                             "%s|%s" % (rule, record["grammar_id"]), kept, allowance)
        found = []
        examined = certified = 0
        for definition in chosen:
            product = B.realise(labelled, definition)
            examined += 1
            if not product.built:
                continue
            cert = B.certify(product.mol, labelled, required, expected, len(ion.core))
            if not cert.certified:
                continue
            certified += 1
            if product.smiles == held_smiles or product.smiles in seen:
                continue
            found.append((definition, product, cert, S.similarity(ion.mol, product.mol),
                          S.similarity(start, product.mol)))
        step = {"round": round_number, "held_smiles": held_smiles, "held_heavy_atoms": int(ion.mol.GetNumAtoms()),
                "candidates_examined": examined, "candidates_certified": certified, "capacity": capacity(kept)}
        if not found:
            step.update({"moved": False, "reason": "no certified product that is new to this chain"})
            steps.append(step)
            break

        def ordering(item):
            key = Q.pick_key(seed, record["family_id"], record["grammar_id"], rule, round_number, item[1].smiles)
            if preference == "most similar":
                return (-item[4], key)
            if preference == "least similar":
                return (item[4], key)
            return (0.0, key)

        found.sort(key=ordering)
        definition, product, cert, near_held, near_start = found[0]
        rec = _record(ion, {**record, "canonical_smiles": held_smiles}, definition, 0, product, cert,
                      ion.mol, query, required, ws, svg=svg, width=width, height=height)
        rec.update({"similarity_to_held": near_held, "similarity_to_start": near_start,
                    "candidates_offered": len(found)})
        step.update({"moved": True, "selected": rec})
        steps.append(step)
        seen.add(product.smiles)
        held = product.mol
        held_core = sorted(int(a.GetIdx()) for a in product.mol.GetAtoms()
                           if a.HasProp(B.CORE) and int(a.GetIntProp(B.CORE)) == 1)

    last = Chem.MolToSmiles(held)
    if progress:
        progress(1.0, "done")
    return {"rule": rule, "preference": preference, "depth": depth, "steps": steps,
            "start": {"smiles": record["canonical_smiles"], "family_id": record["family_id"],
                      "family_label": record.get("family_label", ""), "grammar_id": record["grammar_id"],
                      "core_pattern": record["core_pattern"], "core_atoms": sorted(int(i) for i in str(record["core_atoms"]).split())},
            "summary": {"rounds_moved": sum(1 for s in steps if s.get("moved")), "final_smiles": last,
                        "final_heavy_atoms": int(held.GetNumAtoms()),
                        "similarity_first_to_last": S.similarity(start, held),
                        "elapsed_s": round(time.perf_counter() - t0, 2)},
            "settings": {"seed": int(seed), "allowance": int(allowance), "operators": sorted(ops), "workspace": ws.name}}
