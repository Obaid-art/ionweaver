
# -*- coding: utf-8 -*-
"""Stage 04 — body components, the library, and terminal fragments.

Three things with three names:

  body component   a connected part of one ion's body, as that ion stands
  library component  what a body component becomes once pooled across families,
                   deduplicated, and recorded with its attachment pattern
  terminal fragment  the part of either lying beyond an acyclic bridge bond on
                   the side away from every attachment point

A component is recorded with a dummy atom standing at each attachment, bonded at
the order the attachment carried. The canonical SMILES of that structure is the
library key, so two body components are the same library component only when
their attachment pattern and its bond orders agree as well as their structure.

The body is cut only at acyclic bonds. This follows from the ring principle
and is not a separate rule: cutting a ring bond would break a ring, so a ring
the core covers in part is held whole and its body atoms join the fixed set.
They are not core -- certification does not hold their hydrogen count or their
bonds to anything outside the ring -- and they are not part of any component, so
they are never removed or exchanged as a group. They stay editable wherever no
operator alters a ring, which is an appended insertion or a graft at a free
valence. Without this the half of a fused aromatic ring lying outside the core
would be offered as a component, and a benzimidazolide would be asked to give up
four carbons of its benzo ring.

Eligibility is reported per operator because it differs between them, and a
component excluded from an operator stays in the instances file with the reason.
Insert, Substitute and Rebond do not draw on the library at all.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

from .common import RunError, parse

#: a library component must reach this to be available to the group-level four
MIN_GROUP_HEAVY_ATOMS = 2

#: a terminal fragment must reach this to be prunable or exchangeable
MIN_FRAGMENT_HEAVY_ATOMS = 2

INSTANCE_COLUMNS = (
    "instance_id", "inventory_id", "family_id", "family_label", "core_pattern",
    "component_smiles", "component_key", "heavy_atoms", "atom_indices",
    "attachment_count", "attachment_roots", "attachment_bond_orders",
    "attaches_to_core_atoms", "rings_in_component", "spans_a_ring_with_the_core",
    "admitted", "reason")

LIBRARY_COLUMNS = (
    "library_id", "component_key", "component_smiles", "heavy_atoms",
    "attachment_count", "attachment_bond_orders", "rings_in_component",
    "occurrences", "ions", "families", "first_inventory_id", "first_family_id",
    "terminal_fragments_at_two_or_more",
    "graft_eligible", "graft_reason",
    "prune_eligible", "prune_reason",
    "swap_eligible", "swap_reason",
    "crossover_eligible", "crossover_reason")

FRAGMENT_COLUMNS = (
    "fragment_id", "fragment_key", "fragment_smiles", "heavy_atoms",
    "cut_bond_order", "from_library_id", "from_component_key",
    "is_itself_a_library_component", "occurrences", "source_libraries")


# ---------------------------------------------------------------------------
# cutting


def fixed_atoms(mol: Chem.Mol,
                core_atoms: Sequence[int]) -> Tuple[Set[int], Set[int]]:
    """The atoms no component may contain, and which of them are not core.

    Every ring holding at least one fixed atom is taken whole, repeatedly, so a
    fused system is absorbed as far as the fusion reaches. The second set is
    what the closure added: atoms outside the core that lie in a ring the core
    covers in part.
    """
    keep = set(int(i) for i in core_atoms)
    rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    growing = True
    while growing:
        growing = False
        for ring in rings:
            if ring & keep and not ring <= keep:
                keep |= ring
                growing = True
    return keep, keep - set(int(i) for i in core_atoms)



def _copy_out(mol: Chem.Mol, keep: Sequence[int],
              attachments: Sequence[Tuple[int, Chem.BondType]]) -> Chem.Mol:
    """The atoms of `keep` as their own molecule, with a dummy per attachment.

    `attachments` names, for each bond that left the set, the atom inside it
    came from and the order it carried. A dummy preserves that atom's valence,
    so its hydrogen count is the one it had in the ion.
    """
    work = Chem.RWMol()
    index: Dict[int, int] = {}
    for old in keep:
        index[old] = work.AddAtom(Chem.Atom(mol.GetAtomWithIdx(old)))
    for old in keep:
        for bond in mol.GetAtomWithIdx(old).GetBonds():
            other = bond.GetOtherAtomIdx(old)
            if other in index and other > old:
                new = work.AddBond(index[old], index[other],
                                   bond.GetBondType()) - 1
                work.GetBondWithIdx(new).SetIsAromatic(bond.GetIsAromatic())
    for root, order in attachments:
        dummy = work.AddAtom(Chem.Atom(0))
        work.AddBond(index[root], dummy, order)
    out = work.GetMol()
    Chem.SanitizeMol(out)
    return out


@dataclass
class Instance:
    """One body component of one ion."""
    inventory_id: int
    family_id: str
    atom_indices: Tuple[int, ...]
    attachment_roots: Tuple[int, ...]
    attachment_orders: Tuple[str, ...]
    core_neighbours: Tuple[int, ...]
    spans_a_ring_with_the_core: bool
    smiles: str = ""
    key: str = ""
    rings: int = 0
    admitted: bool = False
    reason: str = ""


def components_of(mol: Chem.Mol, core_atoms: Sequence[int]) -> List[Instance]:
    """Every body component of one ion, with its attachments.

    The body is every atom the fixed set does not hold, and its connected parts
    are the body components. An attachment is a bond from a component atom to a
    fixed atom. Because the fixed set is closed over rings, no attachment is
    ever a ring bond, so every component cuts out as a molecule and no ring is
    ever divided.

    Every component then attaches at exactly one point, and that follows.
    The fixed set is connected: the core is connected, and each ring the closure
    absorbs shares an atom with what is already there. A component attached at
    two points would close a cycle through both attachment bonds and through the
    fixed set, making both of them ring bonds, which the closure has already
    ruled out. So multi-attachment is not a case to exclude; it is a case that
    cannot arise, and the run asserts it below.
    """
    core, _absorbed = fixed_atoms(mol, core_atoms)
    body = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
    if not body:
        return []

    seen: Set[int] = set()
    out: List[Instance] = []
    for start in body:
        if start in seen:
            continue
        group: List[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            here = stack.pop()
            group.append(here)
            for neighbour in mol.GetAtomWithIdx(here).GetNeighbors():
                j = int(neighbour.GetIdx())
                if j not in core and j not in seen:
                    seen.add(j)
                    stack.append(j)
        group.sort()

        attachments: List[Tuple[int, Chem.BondType]] = []
        roots: List[int] = []
        orders: List[str] = []
        neighbours: List[int] = []
        ring_spanning = False
        for inside in group:
            for bond in mol.GetAtomWithIdx(inside).GetBonds():
                other = int(bond.GetOtherAtomIdx(inside))
                if other in core:
                    attachments.append((inside, bond.GetBondType()))
                    roots.append(inside)
                    orders.append(str(bond.GetBondType()))
                    neighbours.append(other)
                    if bond.IsInRing():
                        ring_spanning = True

        instance = Instance(
            inventory_id=0, family_id="",
            atom_indices=tuple(group),
            attachment_roots=tuple(roots),
            attachment_orders=tuple(orders),
            core_neighbours=tuple(neighbours),
            spans_a_ring_with_the_core=ring_spanning)

        if ring_spanning:
            raise RunError(
                "a component of %s is attached by a ring bond, which the ring "
                "closure of the fixed set exists to make impossible: atoms %s"
                % (Chem.MolToSmiles(mol), group))
        if len(roots) != 1:
            raise RunError(
                "a component of %s attaches at %d points. One is the only "
                "possibility once the fixed set is ring-closed and connected, "
                "so this means the closure or the connectivity failed: atoms "
                "%s at %s." % (Chem.MolToSmiles(mol), len(roots), group, roots))
        try:
            piece = _copy_out(mol, group, attachments)
        except (ValueError, Chem.MolSanitizeException) as problem:
            # A component is cut only at acyclic bonds, so a cut breaks nothing.
            # A failure here is a defect in the cutting and not a property of
            # the chemistry, so it stops the run instead of quietly removing a
            # component from the library.
            raise RunError(
                "cutting a component out of %s did not sanitize: %s. Atoms %s, "
                "attached at %s." % (Chem.MolToSmiles(mol),
                                     str(problem).strip().splitlines()[0],
                                     group, roots))

        instance.smiles = Chem.MolToSmiles(piece)
        instance.key = instance.smiles
        instance.rings = piece.GetRingInfo().NumRings()

        charge = sum(int(a.GetFormalCharge()) for a in piece.GetAtoms())
        radicals = sum(int(a.GetNumRadicalElectrons()) for a in piece.GetAtoms())
        if charge:
            instance.reason = ("the component carries a formal charge of %+d; "
                               "all charge belongs to the core" % charge)
        elif radicals:
            instance.reason = ("the component carries %d unpaired electron(s)"
                               % radicals)
        else:
            instance.admitted = True
        out.append(instance)
    return out


# ---------------------------------------------------------------------------
# terminal fragments


def terminal_fragments(piece: Chem.Mol) -> List[Dict[str, Any]]:
    """Every terminal fragment of a component, at two or more heavy atoms.

    A terminal fragment lies beyond an acyclic bridge bond, on the side holding
    no attachment root. The bond cut is always internal to the component: a
    dummy bond is the attachment itself and is never a candidate, which is what
    keeps a cross-over from coinciding with a swap.
    """
    dummies = frozenset(a.GetIdx() for a in piece.GetAtoms()
                        if a.GetAtomicNum() == 0)
    roots = frozenset(
        int(bond.GetOtherAtomIdx(index))
        for index in dummies
        for bond in piece.GetAtomWithIdx(index).GetBonds())
    heavy = frozenset(a.GetIdx() for a in piece.GetAtoms()
                      if a.GetAtomicNum() != 0)

    out: List[Dict[str, Any]] = []
    for cut in terminal_pieces(piece, heavy, roots, dummies):
        fragment = _copy_out(piece, sorted(cut.atoms), [(cut.root, cut.order)])
        out.append({"fragment_key": Chem.MolToSmiles(fragment),
                    "heavy_atoms": len(cut.atoms),
                    "cut_bond_order": str(cut.order)})
    return out


@dataclass(frozen=True)
class Cut:
    """One terminal piece, named from both ends of the bond that frees it.

    `atoms` is the piece. `root` is its own atom at the cut, where a dummy
    stands once the piece is copied out. `anchor` is the atom on the side that
    stays, which is the atom an edit touches. Both are named because they are
    the two ends of one bond and confusing them is easy: the library build needs
    the root and the enumerator needs the anchor.
    """
    atoms: FrozenSet[int]
    root: int
    anchor: int
    order: Chem.BondType


def terminal_pieces(mol: Chem.Mol, atoms: FrozenSet[int],
                    roots: FrozenSet[int],
                    ignore: FrozenSet[int] = frozenset()) -> List[Cut]:
    """Every terminal piece of `atoms`, named from both ends of its cut.

    One entry per acyclic bridge bond with both ends inside `atoms` whose far
    side holds no root, no ignored atom, and at least the fragment floor's worth
    of heavy atoms. The walk never leaves `atoms` together with `ignore`, so a
    structure's attachment bond is never a candidate and a ring is never cut.

    The enumerator calls this on a body component's atoms in the ion, and the
    library build calls it on a component cut out with dummies standing at its
    attachments. Both need the same answer, so there is one function.
    """
    inside = set(atoms) | set(ignore)
    out: List[Cut] = []
    for bond in mol.GetBonds():
        begin, end = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        if begin not in atoms or end not in atoms or bond.IsInRing():
            continue

        sides: List[Set[int]] = []
        for seed, blocked in ((begin, end), (end, begin)):
            reached = {seed}
            stack = [seed]
            while stack:
                here = stack.pop()
                for neighbour in mol.GetAtomWithIdx(here).GetNeighbors():
                    j = int(neighbour.GetIdx())
                    if j not in inside:
                        continue
                    if here in (begin, end) and j == blocked:
                        continue
                    if j not in reached:
                        reached.add(j)
                        stack.append(j)
            sides.append(reached)
        if sides[0] & sides[1]:
            continue  # the bond is in a cycle after all, so it bridges nothing

        # (the piece, its own atom at the cut, the atom that stays)
        for side, root, anchor in ((sides[0], begin, end),
                                   (sides[1], end, begin)):
            if side & roots or side & set(ignore):
                continue
            if len(side) < MIN_FRAGMENT_HEAVY_ATOMS:
                continue
            out.append(Cut(frozenset(side), root, anchor,
                           bond.GetBondType()))
    return out


# ---------------------------------------------------------------------------
# eligibility


def eligibility(heavy_atoms: int, attachment_count: int,
                fragments: int) -> Dict[str, Any]:
    """Which operators may use this library component, and why not where not.

    Graft and Swap are the two that require a single attachment point: a
    component attached at more than one forms a ring with the core, and the
    same attachment count does not imply the same span. Prune and Cross-over
    need a terminal fragment, so a component with none is unavailable to them
    however large it is.
    """
    out: Dict[str, Any] = {}

    small = heavy_atoms < MIN_GROUP_HEAVY_ATOMS
    many = attachment_count != 1
    no_fragment = fragments == 0

    def say(name: str, ok: bool, reason: str) -> None:
        out[name + "_eligible"] = int(ok)
        out[name + "_reason"] = "" if ok else reason

    reason_small = ("one heavy atom; a single-atom addition is an Insert and a "
                    "single-atom exchange a Substitute")
    reason_many = ("attached at %d points; it forms a ring with the core, and "
                   "no operator alters a ring" % attachment_count)
    reason_fragment = ("no terminal fragment of %d or more heavy atoms at an "
                       "acyclic bridge bond" % MIN_FRAGMENT_HEAVY_ATOMS)

    say("graft", not small and not many,
        reason_small if small else reason_many)
    say("swap", not small and not many,
        reason_small if small else reason_many)
    say("prune", not small and not no_fragment,
        reason_small if small else reason_fragment)
    say("crossover", not small and not no_fragment,
        reason_small if small else reason_fragment)
    return out


# ---------------------------------------------------------------------------
# the stage


def build(starting_ions: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """The instances, the library and the terminal fragments.

    `starting_ions` carries, per editable ion, its canonical SMILES, its family
    and its core atom indices.
    """
    instances: List[Dict[str, Any]] = []
    pooled: Dict[str, Dict[str, Any]] = {}
    held: List[Dict[str, Any]] = []
    instance_id = 0

    for record in starting_ions:
        mol = parse(record["canonical_smiles"])
        core = [int(i) for i in str(record["core_atoms"]).split()]
        _fixed, absorbed = fixed_atoms(mol, core)
        if absorbed:
            held.append({
                "inventory_id": record["inventory_id"],
                "family_id": record["family_id"],
                "canonical_smiles": record["canonical_smiles"],
                "core_pattern": record.get("core_pattern", ""),
                "core_atoms": " ".join(str(i) for i in sorted(core)),
                "ring_shared_atoms": " ".join(str(i) for i
                                              in sorted(absorbed)),
                "ring_shared_count": len(absorbed),
                "heavy_atoms": mol.GetNumAtoms(),
                "reason": "these atoms lie in a ring the core covers in part, "
                          "so the ring is held whole; they are not core and "
                          "they are in no component, and an appended insertion "
                          "or a graft at a free valence still reaches them"})
        for found in components_of(mol, core):
            instance_id += 1
            instances.append({
                "instance_id": instance_id,
                "inventory_id": record["inventory_id"],
                "family_id": record["family_id"],
                "family_label": record.get("family_label", ""),
                "core_pattern": record.get("core_pattern", ""),
                "component_smiles": found.smiles,
                "component_key": found.key,
                "heavy_atoms": len(found.atom_indices),
                "atom_indices": " ".join(str(i) for i in found.atom_indices),
                "attachment_count": len(found.attachment_roots),
                "attachment_roots": " ".join(str(i) for i
                                             in found.attachment_roots),
                "attachment_bond_orders": " ".join(found.attachment_orders),
                "attaches_to_core_atoms": " ".join(str(i) for i
                                                   in found.core_neighbours),
                "rings_in_component": found.rings,
                "spans_a_ring_with_the_core":
                    int(found.spans_a_ring_with_the_core),
                "admitted": int(found.admitted),
                "reason": found.reason})

            if not found.admitted:
                continue
            entry = pooled.setdefault(found.key, {
                "component_key": found.key,
                "component_smiles": found.smiles,
                "heavy_atoms": len(found.atom_indices),
                "attachment_count": len(found.attachment_roots),
                "attachment_bond_orders": " ".join(found.attachment_orders),
                "rings_in_component": found.rings,
                "occurrences": 0,
                "ions": set(), "families": set(),
                "first_inventory_id": record["inventory_id"],
                "first_family_id": record["family_id"]})
            entry["occurrences"] += 1
            entry["ions"].add(record["inventory_id"])
            entry["families"].add(record["family_id"])

    # the library, in a stable order: size, then attachments, then the key
    library: List[Dict[str, Any]] = []
    fragment_rows: Dict[str, Dict[str, Any]] = {}
    order = sorted(pooled.values(),
                   key=lambda e: (e["heavy_atoms"], e["attachment_count"],
                                  e["component_key"]))
    for number, entry in enumerate(order, start=1):
        piece = parse(entry["component_key"])
        found = terminal_fragments(piece)
        row = dict(entry)
        row["library_id"] = "L%04d" % number
        row["ions"] = len(entry["ions"])
        row["families"] = len(entry["families"])
        row["terminal_fragments_at_two_or_more"] = len(found)
        row.update(eligibility(entry["heavy_atoms"], entry["attachment_count"],
                               len(found)))
        library.append(row)

        for fragment in found:
            seen = fragment_rows.setdefault(fragment["fragment_key"], {
                "fragment_key": fragment["fragment_key"],
                "fragment_smiles": fragment["fragment_key"],
                "heavy_atoms": fragment["heavy_atoms"],
                "cut_bond_order": fragment["cut_bond_order"],
                "from_library_id": row["library_id"],
                "from_component_key": entry["component_key"],
                "occurrences": 0, "sources": set()})
            seen["occurrences"] += 1
            seen["sources"].add(row["library_id"])

    keys = {row["component_key"] for row in library}
    fragments: List[Dict[str, Any]] = []
    for number, fragment in enumerate(
            sorted(fragment_rows.values(),
                   key=lambda f: (f["heavy_atoms"], f["fragment_key"])),
            start=1):
        sources = fragment.pop("sources")
        fragment["fragment_id"] = "T%04d" % number
        fragment["source_libraries"] = " ".join(sorted(sources))
        fragment["is_itself_a_library_component"] = int(
            fragment["fragment_key"] in keys)
        fragments.append(fragment)

    for row in library:
        row.pop("component_smiles_set", None)

    return {"instances": instances, "library": library, "fragments": fragments,
            "ring_shared": held,
            "settings": {"min_group_heavy_atoms": MIN_GROUP_HEAVY_ATOMS,
                         "min_fragment_heavy_atoms": MIN_FRAGMENT_HEAVY_ATOMS,
                         "library_key": "canonical SMILES of the component "
                                        "with a dummy atom at each attachment, "
                                        "bonded at the order the attachment "
                                        "carried",
                         "body_cut_only_at_acyclic_bonds": True}}


def summary(result: Dict[str, Any]) -> Dict[str, Any]:
    """Counts, every one computed from the rows."""
    instances = result["instances"]
    library = result["library"]
    fragments = result["fragments"]
    not_admitted = collections.Counter(
        row["reason"] for row in instances if not row["admitted"])
    return {
        "component_instances": len(instances),
        "instances_admitted": sum(row["admitted"] for row in instances),
        "instances_not_admitted": dict(not_admitted),
        "library_components": len(library),
        "single_atom_components": sum(1 for row in library
                                      if row["heavy_atoms"] == 1),
        "multi_attachment_components": sum(1 for row in library
                                           if row["attachment_count"] > 1),
        "multi_attachment_occurrences": sum(
            1 for row in instances
            if row["admitted"] and row["attachment_count"] > 1),
        "ions_with_a_multi_attachment_component": len({
            row["inventory_id"] for row in instances
            if row["admitted"] and row["attachment_count"] > 1}),
        "components_attached_by_a_ring_bond": sum(
            row["spans_a_ring_with_the_core"] for row in instances),
        "ions_with_ring_shared_atoms": len(result["ring_shared"]),
        "ring_shared_atoms": sum(row["ring_shared_count"]
                                 for row in result["ring_shared"]),
        "largest_component_heavy_atoms": max(
            [row["heavy_atoms"] for row in library] or [0]),
        "eligible": {name: sum(row[name + "_eligible"] for row in library)
                     for name in ("graft", "prune", "swap", "crossover")},
        "terminal_fragments_at_two_or_more": len(fragments),
        "fragments_that_are_not_library_components": sum(
            1 for row in fragments
            if not row["is_itself_a_library_component"]),
    }


def starting_ions(cores: Dict[str, Any],
                  inventory: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """The editable set: every ion a core reached, with its fixed atoms.

    One row per assignment, carrying the core atoms, the atoms the ring
    closure then holds, and the size of the body that is
    left. An ion whose body is empty has a core and no editable position; it is
    reported as such and never as an ion without a core.
    """
    records = {r["inventory_id"]: r for r in inventory}
    pattern = {m["family_id"] + "|" + m["grammar_id"]: m["pattern"]
               for m in cores["models"]}
    out: List[Dict[str, Any]] = []
    for row in cores["assignments"]:
        record = records[row["inventory_id"]]
        mol = parse(record["canonical_smiles"])
        core = sorted(int(i) for i in row["matched_atoms"].split())
        fixed, absorbed = fixed_atoms(mol, core)
        pieces = components_of(mol, core)
        out.append({
            "inventory_id": row["inventory_id"],
            "family_id": row["family_id"],
            "family_label": record["family_label"],
            "canonical_smiles": record["canonical_smiles"],
            "grammar_id": row["grammar_id"],
            "core_pattern": pattern[row["family_id"] + "|" + row["grammar_id"]],
            "core_atoms": " ".join(str(i) for i in core),
            "core_atom_count": len(core),
            "ring_shared_atoms": " ".join(str(i) for i in sorted(absorbed)),
            "ring_shared_count": len(absorbed),
            "heavy_atoms": record["heavy_atoms"],
            "body_heavy_atoms": mol.GetNumAtoms() - len(fixed),
            "body_components": len(pieces)})
    return out
