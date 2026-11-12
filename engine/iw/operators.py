
# -*- coding: utf-8 -*-
"""Stage 05a — the eight operators, enumerated as definitions.

Spec B and C. A definition is a specified edit before anything is built: what it
removes, what it puts back, and where. Nothing here assembles a molecule. The
conditions of C1 are applied to definitions, so an edit that cannot be given a
legal valence state, or whose touched atoms would land outside the sets the
inventory displays, is never created and consumes nothing.

The signature is what makes B0's exclusion principle testable. Every
definition is keyed by the change it makes to the parent graph, in a form that
carries no operator's name: the atoms it removes, the bonds it removes, the
bonds it adds between atoms already there, the elements it re-specifies, the
orders it re-specifies, and the structures it attaches with the atom and order
each attaches at. Two operators reaching the same signature on the same parent
would be the same edit under two names, and that stops the run.

A signature names the parent's own atom indices, so two definitions on different
ions can never collide. The test therefore runs per ion and the signatures are
discarded after it, which is what makes an exhaustive test over several million
definitions possible at all.

One clarification the specification left open. C1 governs an atom the edit
*constructs* -- Insert's new atom and Substitute's new element -- which must land
in a creatable environment: observed, neutral and non-aromatic. An atom arriving
inside a fragment cut from a curated ion is not constructed; it is carried, and
it arrives in the environment it already occupied there. Read the other way, no
aromatic component could ever be grafted, which would contradict B3's whole
account of how rings enter a structure. Carried atoms are asserted to be in
observed environments, not filtered against the creatable set.
"""
from __future__ import annotations

import collections
from dataclasses import dataclass, field
from typing import (Any, Dict, FrozenSet, Iterable, Iterator, List, Optional,
                    Sequence, Set, Tuple)

from rdkit import Chem

from .common import (AtomEnv, BondEnv, RunError, atom_environment,
                     bond_end, bond_key, parse)
from .fragments import MIN_FRAGMENT_HEAVY_ATOMS, MIN_GROUP_HEAVY_ATOMS, Cut, \
    terminal_pieces
from .ledger import BOND_TYPE, ORDER_VALUE, Placed, Resulting, Touch, placed, \
    retained
from .registry import RESTRICTIONS, REASONS, restriction_failed

PT = Chem.GetPeriodicTable()

#: The orders an edit may put a bond at. No aromatic order appears, because no
#: operator alters a ring bond.
ORDERS: Tuple[str, ...] = ("SINGLE", "DOUBLE", "TRIPLE")

ATOM_LEVEL = ("Insert", "Substitute", "Rebond", "Delete")
GROUP_LEVEL = ("Graft", "Prune", "Swap", "Cross-over")
#: Every operator and the modes it is reported under. A mode is a case the
#: operator separates for reporting; it is never a second operator.
MODES: Dict[str, Tuple[str, ...]] = {
    "Insert": ("bridged", "appended"),
    "Substitute": ("element",),
    "Rebond": ("order",),
    "Delete": ("degree_two", "degree_one"),
    "Graft": ("attach",),
    "Prune": ("cut",),
    "Swap": ("exchange",),
    "Cross-over": ("splice",),
}

#: The stages a definition passes, in the order C4 reports them.
CONDITIONS: Tuple[str, ...] = (("ledger", "atom_environment",
                                "bond_environment") + RESTRICTIONS)


# ---------------------------------------------------------------------------
# describing a structure an edit attaches


_KEY_CACHE: Dict[Tuple[int, Tuple[str, ...]], str] = {}


def built_atom_key(atomic_number: int, orders: Sequence[str]) -> str:
    """One constructed atom, written the way a library component is written.

    A dummy stands at each bond the atom will make, at that bond's order, so an
    atom Insert builds and a component Graft attaches are described in the same
    language. That is what lets one signature compare them, and it is how the
    floor keeping single-atom additions out of Graft is tested rather than
    assumed.

    This describes; it does not validate. A fluorine bridged into a bond
    would carry two bonds, which no valence allows, and that definition still
    has to be offered so the ledger can be the thing that declines it and the
    shrinkage table can say so. Writing the key is therefore done without
    sanitization, which would raise on exactly the cases the method is counting.
    """
    cached = _KEY_CACHE.get((atomic_number, tuple(orders)))
    if cached is not None:
        return cached
    work = Chem.RWMol()
    centre = work.AddAtom(Chem.Atom(int(atomic_number)))
    for order in orders:
        dummy = work.AddAtom(Chem.Atom(0))
        work.AddBond(centre, dummy, BOND_TYPE[order])
    mol = work.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    key = Chem.MolToSmiles(mol)
    _KEY_CACHE[(atomic_number, tuple(orders))] = key
    return key


# ---------------------------------------------------------------------------
# a definition


@dataclass(frozen=True)
class Attachment:
    """A structure an edit joins on, and where it joins."""
    structure: str
    anchors: Tuple[Tuple[int, str], ...]


@dataclass
class Definition:
    """One specified edit, before anything is built."""
    operator: str
    mode: str
    remove_atoms: FrozenSet[int] = frozenset()
    remove_bonds: FrozenSet[Tuple[int, int]] = frozenset()
    add_bonds: FrozenSet[Tuple[int, int, str]] = frozenset()
    element_changes: Tuple[Tuple[int, int], ...] = ()
    order_changes: Tuple[Tuple[int, int, str], ...] = ()
    attachments: Tuple[Attachment, ...] = ()
    touches: Tuple[Touch, ...] = ()
    builds: Tuple[Placed, ...] = ()
    site: str = ""
    description: str = ""
    heavy_atom_change: int = 0
    #: filled by the conditions
    condition_failed: str = ""
    condition_reason: str = ""
    condition_atom: int = -1
    results: Tuple[Resulting, ...] = ()

    @property
    def signature(self) -> Tuple[Any, ...]:
        """The change this edit makes, with no operator's name in it."""
        return (
            tuple(sorted(self.remove_atoms)),
            tuple(sorted(self.remove_bonds)),
            tuple(sorted(self.add_bonds)),
            tuple(sorted(self.element_changes)),
            tuple(sorted(self.order_changes)),
            tuple(sorted((a.structure, a.anchors) for a in self.attachments)),
        )

    @property
    def signature_text(self) -> str:
        atoms, bonds, added, elements, orders, attached = self.signature
        parts = []
        if atoms:
            parts.append("remove atoms " + " ".join(str(i) for i in atoms))
        if bonds:
            parts.append("remove bonds "
                         + " ".join("%d-%d" % b for b in bonds))
        if added:
            parts.append("add bonds "
                         + " ".join("%d-%d:%s" % b for b in added))
        if elements:
            parts.append("set element "
                         + " ".join("%d:%s" % (i, PT.GetElementSymbol(z))
                                    for i, z in elements))
        if orders:
            parts.append("set order "
                         + " ".join("%d-%d:%s" % o for o in orders))
        for structure, anchors in attached:
            parts.append("attach %s at %s"
                         % (structure,
                            " ".join("%d:%s" % a for a in anchors)))
        return "; ".join(parts)


# ---------------------------------------------------------------------------
# one ion's sites


@dataclass
class Component:
    """One body component of the ion being edited."""
    atoms: FrozenSet[int]
    root: int
    anchor: int
    order: str
    key: str


@dataclass
class Context:
    """Everything one starting ion offers, resolved once."""
    inventory_id: int
    family_id: str
    canonical_smiles: str
    mol: Chem.Mol
    core: FrozenSet[int]
    fixed: FrozenSet[int]
    components: Tuple[Component, ...]
    elements: Tuple[int, ...]

    def bond_key(self, a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a < b else (b, a)

    def editable_bond(self, bond: Chem.Bond) -> bool:
        """An acyclic bond that is not internal to the core.

        Internal means both ends in the core. A bond from the core out to the
        body is editable: it is not part of the core's induced bond set, so
        certification does not hold it, and changing it leaves the core atom's
        element, charge and aromaticity alone.
        """
        if bond.IsInRing():
            return False
        a, b = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        return not (a in self.core and b in self.core)

    def body_atom(self, index: int) -> bool:
        return index not in self.fixed


def context(record: Dict[str, Any], elements: Sequence[int]) -> Context:
    """Resolve one starting ion's sites from its row in `starting_ions`."""
    from .fragments import components_of, fixed_atoms
    mol = parse(record["canonical_smiles"])
    core = frozenset(int(i) for i in str(record["core_atoms"]).split())
    fixed, _absorbed = fixed_atoms(mol, core)
    pieces: List[Component] = []
    for found in components_of(mol, sorted(core)):
        root = found.attachment_roots[0]
        anchor = found.core_neighbours[0]
        pieces.append(Component(frozenset(found.atom_indices), root, anchor,
                                found.attachment_orders[0], found.key))
    return Context(
        inventory_id=record["inventory_id"], family_id=record["family_id"],
        canonical_smiles=record["canonical_smiles"], mol=mol, core=core,
        fixed=frozenset(fixed), components=tuple(pieces),
        elements=tuple(elements))


# ---------------------------------------------------------------------------
# the eight enumerators
#
# Each yields definitions without applying any condition. The conditions run
# once, in one place, so no operator can apply a different set of them.


def insert_bridged(ion: Context) -> Iterator[Definition]:
    """One atom into an acyclic bond, remade at the same order on both sides."""
    for bond in ion.mol.GetBonds():
        if not ion.editable_bond(bond):
            continue
        a, b = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        order = str(bond.GetBondType())
        if order not in ORDER_VALUE:
            continue
        for element in ion.elements:
            yield Definition(
                operator="Insert", mode="bridged",
                remove_bonds=frozenset({ion.bond_key(a, b)}),
                attachments=(Attachment(
                    built_atom_key(element, (order, order)),
                    tuple(sorted(((a, order), (b, order))))),),
                touches=(Touch(a, removed=(order,), added=(order,)),
                         Touch(b, removed=(order,), added=(order,))),
                builds=(Placed(element, (order, order)),),
                site="bond %d-%d" % (a, b),
                description="insert %s into the %s bond between atoms %d and "
                            "%d, remade at that order on both sides"
                            % (PT.GetElementSymbol(element), order.lower(),
                               a, b),
                heavy_atom_change=1)


def insert_appended(ion: Context) -> Iterator[Definition]:
    """One atom onto an existing atom, at order 1, 2 or 3."""
    for atom in ion.mol.GetAtoms():
        index = int(atom.GetIdx())
        for order in ORDERS:
            for element in ion.elements:
                yield Definition(
                    operator="Insert", mode="appended",
                    attachments=(Attachment(
                        built_atom_key(element, (order,)),
                        ((index, order),)),),
                    touches=(Touch(index, added=(order,)),),
                    builds=(Placed(element, (order,)),),
                    site="atom %d" % index,
                    description="append %s to atom %d at a %s bond"
                                % (PT.GetElementSymbol(element), index,
                                   order.lower()),
                    heavy_atom_change=1)


def substitute(ion: Context) -> Iterator[Definition]:
    """One body atom's element, with every bond and order kept.

    A ring atom is excluded: an element change inside a ring can aromatise it
    with no bond touched, and a carbon exchanged for a sulfur in a five-ring
    carrying two double bonds gives thiophene.

    Only the atom itself is touched. An atom environment is that atom's element,
    charge, aromatic status, total valence, attached hydrogens and bond orders;
    it does not record what a neighbour is made of. An element substitution
    keeps every bond and every order, so no neighbour's environment changes and
    there is nothing at a neighbour for a condition to test.
    """
    for atom in ion.mol.GetAtoms():
        index = int(atom.GetIdx())
        if not ion.body_atom(index) or atom.IsInRing():
            continue
        if int(atom.GetFormalCharge()) != 0 or atom.GetIsAromatic():
            continue
        current = int(atom.GetAtomicNum())
        for element in ion.elements:
            if element == current:
                continue
            yield Definition(
                operator="Substitute", mode="element",
                element_changes=((index, element),),
                touches=(Touch(index, new_element=element),),
                site="atom %d" % index,
                description="change atom %d from %s to %s, keeping every bond"
                            % (index, PT.GetElementSymbol(current),
                               PT.GetElementSymbol(element)),
                heavy_atom_change=0)


def rebond(ion: Context) -> Iterator[Definition]:
    """One acyclic bond's order, with both endpoints' elements kept."""
    for bond in ion.mol.GetBonds():
        if not ion.editable_bond(bond):
            continue
        a, b = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        current = str(bond.GetBondType())
        if current not in ORDER_VALUE:
            continue
        for order in ORDERS:
            if order == current:
                continue
            key = ion.bond_key(a, b)
            yield Definition(
                operator="Rebond", mode="order",
                order_changes=((key[0], key[1], order),),
                touches=(Touch(a, removed=(current,), added=(order,)),
                         Touch(b, removed=(current,), added=(order,))),
                site="bond %d-%d" % (a, b),
                description="change the bond between atoms %d and %d from %s "
                            "to %s" % (a, b, current.lower(), order.lower()),
                heavy_atom_change=0)


def delete(ion: Context) -> Iterator[Definition]:
    """One neutral, non-aromatic body atom.

    Two cases, separated by the atom's degree. A degree-two atom not in a ring
    has its neighbours bonded at their shared order; where the two bonds differ
    in order there is no shared order and the case is not defined. A degree-one
    atom is removed and nothing is rejoined.
    """
    for atom in ion.mol.GetAtoms():
        index = int(atom.GetIdx())
        if not ion.body_atom(index):
            continue
        if int(atom.GetFormalCharge()) != 0 or atom.GetIsAromatic():
            continue
        bonds = list(atom.GetBonds())

        if len(bonds) == 1 and not atom.IsInRing():
            bond = bonds[0]
            other = int(bond.GetOtherAtomIdx(index))
            order = str(bond.GetBondType())
            if order not in ORDER_VALUE:
                continue
            yield Definition(
                operator="Delete", mode="degree_one",
                remove_atoms=frozenset({index}),
                remove_bonds=frozenset({ion.bond_key(index, other)}),
                touches=(Touch(other, removed=(order,)),),
                site="atom %d" % index,
                description="remove atom %d, which carries one bond, and "
                            "rejoin nothing" % index,
                heavy_atom_change=-1)

        elif len(bonds) == 2 and not atom.IsInRing():
            first, second = bonds
            x = int(first.GetOtherAtomIdx(index))
            y = int(second.GetOtherAtomIdx(index))
            ox, oy = str(first.GetBondType()), str(second.GetBondType())
            if ox != oy or ox not in ORDER_VALUE:
                continue
            if x in ion.core and y in ion.core:
                # Unreachable: the fixed set is connected, so two fixed
                # neighbours would put this atom in a ring and the ring
                # exclusion above would already have skipped it. Kept as the
                # assertion that says so.
                raise RunError(
                    "atom %d of %s is not in a ring and has two core "
                    "neighbours, which the connected fixed set makes "
                    "impossible" % (index, ion.canonical_smiles))
            key = ion.bond_key(x, y)
            yield Definition(
                operator="Delete", mode="degree_two",
                remove_atoms=frozenset({index}),
                remove_bonds=frozenset({ion.bond_key(index, x),
                                        ion.bond_key(index, y)}),
                add_bonds=frozenset({(key[0], key[1], ox)}),
                touches=(Touch(x, removed=(ox,), added=(ox,)),
                         Touch(y, removed=(oy,), added=(oy,))),
                site="atom %d" % index,
                description="remove atom %d and bond its neighbours %d and %d "
                            "at their shared %s order"
                            % (index, x, y, ox.lower()),
                heavy_atom_change=-1)


def graft(ion: Context, library: Sequence[Dict[str, Any]]
          ) -> Iterator[Definition]:
    """A library component of two or more heavy atoms, at a position that
    carried none. The anchor may be a core atom."""
    for atom in ion.mol.GetAtoms():
        index = int(atom.GetIdx())
        for entry in library:
            order = entry["attachment_bond_orders"]
            yield Definition(
                operator="Graft", mode="attach",
                attachments=(Attachment(entry["component_key"],
                                        ((index, order),)),),
                touches=(Touch(index, added=(order,)),),
                site="atom %d" % index,
                description="attach %s to atom %d at a %s bond"
                            % (entry["component_key"], index, order.lower()),
                heavy_atom_change=entry["heavy_atoms"])


def prune(ion: Context) -> Iterator[Definition]:
    """A terminal fragment of two or more heavy atoms, at an acyclic bridge
    bond inside a body component, holding no attachment root."""
    for piece in ion.components:
        for cut in terminal_pieces(ion.mol, piece.atoms,
                                   frozenset({piece.root})):
            order = str(cut.order)
            if order not in ORDER_VALUE:
                continue
            yield Definition(
                operator="Prune", mode="cut",
                remove_atoms=cut.atoms,
                remove_bonds=frozenset({ion.bond_key(cut.anchor, cut.root)}),
                touches=(Touch(cut.anchor, removed=(order,)),),
                site="bond %d-%d" % (cut.anchor, cut.root),
                description="remove the %d-atom fragment beyond the %s bond "
                            "between atoms %d and %d"
                            % (len(cut.atoms), order.lower(), cut.anchor,
                               cut.root),
                heavy_atom_change=-len(cut.atoms))


def swap(ion: Context, library: Sequence[Dict[str, Any]]
         ) -> Iterator[Definition]:
    """A whole body component of two or more heavy atoms for a different
    library component presenting the same bond order at its attachment."""
    for piece in ion.components:
        if len(piece.atoms) < MIN_GROUP_HEAVY_ATOMS:
            continue
        for entry in library:
            if entry["attachment_bond_orders"] != piece.order:
                continue
            if entry["component_key"] == piece.key:
                continue
            yield Definition(
                operator="Swap", mode="exchange",
                remove_atoms=piece.atoms,
                remove_bonds=frozenset({ion.bond_key(piece.anchor,
                                                     piece.root)}),
                attachments=(Attachment(entry["component_key"],
                                        ((piece.anchor, piece.order),)),),
                touches=(Touch(piece.anchor, removed=(piece.order,),
                               added=(piece.order,)),),
                site="component at atom %d" % piece.anchor,
                description="exchange the %d-atom component %s at atom %d for "
                            "%s" % (len(piece.atoms), piece.key, piece.anchor,
                                    entry["component_key"]),
                heavy_atom_change=entry["heavy_atoms"] - len(piece.atoms))


def crossover(ion: Context, fragments: Sequence[Dict[str, Any]]
              ) -> Iterator[Definition]:
    """A terminal fragment of a body component for a different terminal
    fragment of a library component, at the same bond order.

    The cut bond is never the component's attachment bond to the core, so the
    component always keeps at least one heavy atom beyond that attachment and a
    cross-over can never coincide with a swap.
    """
    from .fragments import _copy_out
    for piece in ion.components:
        for cut in terminal_pieces(ion.mol, piece.atoms,
                                   frozenset({piece.root})):
            order = str(cut.order)
            if order not in ORDER_VALUE:
                continue
            removed = Chem.MolToSmiles(
                _copy_out(ion.mol, sorted(cut.atoms), [(cut.root, cut.order)]))
            for entry in fragments:
                if entry["cut_bond_order"] != order:
                    continue
                if entry["fragment_key"] == removed:
                    continue
                yield Definition(
                    operator="Cross-over", mode="splice",
                    remove_atoms=cut.atoms,
                    remove_bonds=frozenset({ion.bond_key(cut.anchor,
                                                         cut.root)}),
                    attachments=(Attachment(entry["fragment_key"],
                                            ((cut.anchor, order),)),),
                    touches=(Touch(cut.anchor, removed=(order,),
                                   added=(order,)),),
                    site="bond %d-%d" % (cut.anchor, cut.root),
                    description="cut the %d-atom fragment %s beyond the bond "
                                "between atoms %d and %d and join %s in its "
                                "place" % (len(cut.atoms), removed, cut.anchor,
                                           cut.root, entry["fragment_key"]),
                    heavy_atom_change=(entry["heavy_atoms"]
                                       - len(cut.atoms)))




#: the end an attached structure presents, cached by the key that names it
_ATTACHED_END: Dict[str, Any] = {}


def attached_end(structure: str):
    """What a structure joined on presents at its point of attachment.

    A component is written with a dummy atom where it joins, so the end the
    parent will see is the dummy's one neighbour.
    """
    if structure not in _ATTACHED_END:
        mol = Chem.MolFromSmiles(structure, sanitize=False)
        end = None
        if mol is not None:
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetDegree() == 1:
                    end = bond_end(atom.GetNeighbors()[0])
                    break
        _ATTACHED_END[structure] = end
    return _ATTACHED_END[structure]


def bonds_after(ion: Context, definition: Definition,
                results: Sequence[Resulting]) -> List[Tuple[BondEnv, int]]:
    """Every bond the edit makes new, with the atom to blame if it is refused.

    A bond is new when the edit creates it, changes its order, or changes what
    one of its ends is. A bond both of whose ends the edit leaves alone is the
    bond it already was, and was already observed.

    The bonds inside a component joined on are not listed: the component was
    taken whole out of a curated ion, so every bond within it is one that ion
    displayed.
    """
    after: Dict[int, Any] = {}
    for result in results:
        z, charge, aromatic = result.environment[0:3]
        after[result.index] = (int(z), int(charge), bool(aromatic))

    def end_of(index: int):
        if index in after:
            return after[index]
        if index < 0 or index >= ion.mol.GetNumAtoms():
            return None
        return bond_end(ion.mol.GetAtomWithIdx(int(index)))

    out: List[Tuple[BondEnv, int]] = []

    # a structure joined on, at each atom it anchors to
    built = [r for r in results if r.index < 0]
    for attachment in definition.attachments:
        far = attached_end(attachment.structure)
        if far is None and built:
            far = after.get(built[0].index)
        for anchor, order in attachment.anchors:
            near = end_of(int(anchor))
            if near is None or far is None:
                continue
            out.append((bond_key(near, far, order), int(anchor)))

    # a bond the edit draws between two atoms that are already there
    for first, second, order in definition.add_bonds:
        near, far = end_of(int(first)), end_of(int(second))
        if near and far:
            out.append((bond_key(near, far, order), int(first)))

    # a bond whose order the edit changes
    for first, second, order in definition.order_changes:
        near, far = end_of(int(first)), end_of(int(second))
        if near and far:
            out.append((bond_key(near, far, order), int(first)))

    # an atom whose element the edit changes keeps its bonds, but every one of
    # them is now a different pair
    gone = {int(i) for i in definition.remove_atoms}
    cut = {tuple(sorted(b)) for b in definition.remove_bonds}
    for index, _element in definition.element_changes:
        atom = ion.mol.GetAtomWithIdx(int(index))
        for bond in atom.GetBonds():
            other = int(bond.GetOtherAtomIdx(int(index)))
            if other in gone or tuple(sorted((int(index), other))) in cut:
                continue
            near, far = end_of(int(index)), end_of(other)
            if near and far:
                out.append((bond_key(near, far, str(bond.GetBondType())),
                            int(index)))
    return out


# ---------------------------------------------------------------------------
# the conditions, applied once, in one place


def apply_conditions(ion: Context, definition: Definition,
                     observed: FrozenSet[AtomEnv],
                     observed_bonds: FrozenSet[BondEnv]) -> Definition:
    """Decide one definition from the parent graph and the change it specifies.

    The order is the order C4 reports the shrinkage in, so the condition
    recorded against a declined definition is the first one that declined it
    and the columns of the shrinkage table nest.

    An atom the edit constructs answers to the creatable set, which is the
    observed set with the restrictions of C1 applied. An atom it keeps and
    touches answers to the observed set alone, because a core atom a graft
    lands on is aromatic and may be charged, and it already was.
    """
    results: List[Resulting] = []

    for touch in definition.touches:
        results.append(retained(ion.mol, touch))
    for index, build in enumerate(definition.builds):
        results.append(placed(build, -1 - index))
    definition.results = tuple(results)

    for result in results:
        if not result.legal:
            definition.condition_failed = "ledger"
            definition.condition_reason = result.reason
            definition.condition_atom = result.index
            return definition

    for result in results:
        if result.environment not in observed:
            definition.condition_failed = "atom_environment"
            definition.condition_reason = REASONS["not_observed"]
            definition.condition_atom = result.index
            return definition

    if not observed_bonds:
        raise RunError("the bond-environment condition was given no observed "
                       "pairs. An empty set would silently skip the "
                       "condition and quietly run a different method.")
    for env, blame in bonds_after(ion, definition, results):
            if env not in observed_bonds:
                definition.condition_failed = "bond_environment"
                definition.condition_reason = REASONS["bond_not_observed"]
                definition.condition_atom = blame
                return definition

    for result in results:
        if not result.constructed:
            continue
        code, reason = restriction_failed(result.environment)
        if code:
            definition.condition_failed = code
            definition.condition_reason = reason
            definition.condition_atom = result.index
            return definition
    return definition


# ---------------------------------------------------------------------------
# one ion, every operator


def definitions_for(ion: Context, library: Sequence[Dict[str, Any]],
                    fragments: Sequence[Dict[str, Any]],
                    observed: FrozenSet[AtomEnv],
                    observed_bonds: FrozenSet[BondEnv]
                    ) -> Tuple[List[Definition],
                               Dict[Tuple[str, str], Dict[str, int]],
                               List[Dict[str, Any]]]:
    """Every definition the eight operators offer on one ion, decided.

    Returns the definitions that survived every condition, the per-mode tally
    of how many were available at each stage, which is the row C4's shrinkage
    table is built from, and the refusals grouped by what was refused.

    The refusals are grouped because there are over a million of them across
    the run and they differ only in which site they sat at. One row per
    (operator, mode, condition, the environment that was not reachable) with a
    count and one example keeps every distinct piece of chemistry visible
    without listing the same refusal a thousand times.
    """
    graft_library = [e for e in library if e["graft_eligible"]]
    swap_library = [e for e in library if e["swap_eligible"]]

    # Every stream, in one list, so the conditions below are applied to all of
    # them by the same code. Each definition names its own operator and mode.
    sources: List[Iterable[Definition]] = [
        insert_bridged(ion), insert_appended(ion), substitute(ion),
        rebond(ion), delete(ion), graft(ion, graft_library), prune(ion),
        swap(ion, swap_library), crossover(ion, fragments),
    ]

    tally: Dict[Tuple[str, str], Dict[str, int]] = {}
    kept: List[Definition] = []

    def blank() -> Dict[str, int]:
        # `offered` counts every site-and-option pair before any condition, so
        # a definition the ledger refuses outright is still visible. C4's first
        # column is `surviving_after_ledger`, which is the number available
        # under the valence ledger alone.
        row = {"offered": 0}
        row.update({name: 0 for name in CONDITIONS})
        row["accepted"] = 0
        return row

    groups: Dict[Tuple[Any, ...], Dict[str, Any]] = {}

    for stream in sources:
        for definition in stream:
            row = tally.setdefault((definition.operator, definition.mode),
                                   blank())
            row["offered"] += 1
            decided = apply_conditions(ion, definition, observed,
                                       observed_bonds)
            if decided.condition_failed:
                row[decided.condition_failed] += 1
                blame = next((r for r in decided.results
                              if r.index == decided.condition_atom),
                             decided.results[0] if decided.results else None)
                environment = ""
                if blame is not None:
                    z, charge, aromatic, valence, hydrogens, orders =                         blame.environment
                    environment = ("%s charge %+d%s valence %d %dH bonds %s"
                                   % (PT.GetElementSymbol(z), charge,
                                      " aromatic" if aromatic else "",
                                      valence, hydrogens, " ".join(orders)))
                key = (definition.operator, definition.mode,
                       decided.condition_failed, decided.condition_reason,
                       environment)
                found = groups.get(key)
                if found is None:
                    groups[key] = {
                        "inventory_id": ion.inventory_id,
                        "family_id": ion.family_id,
                        "operator": definition.operator,
                        "mode": definition.mode,
                        "condition": decided.condition_failed,
                        "reason": decided.condition_reason,
                        "resulting_environment": environment,
                        "definitions_refused": 1,
                        "example_site": definition.site,
                        "example_description": definition.description}
                else:
                    found["definitions_refused"] += 1
            else:
                row["accepted"] += 1
                kept.append(decided)

    # The surviving count after each condition, which is what C4 asks for: one
    # column per condition, so a reader sees which one does the work.
    for row in tally.values():
        remaining = row["offered"]
        for name in CONDITIONS:
            remaining -= row[name]
            row["surviving_after_" + name] = remaining

    return kept, tally, [groups[k] for k in sorted(groups)]


def disjointness(kept: Sequence[Definition], ion: Context) -> List[Dict[str, Any]]:
    """Assert no signature appears under two operators on this ion.

    A signature names this ion's own atom indices, so a collision can only ever
    be within one ion. Testing per ion and discarding is what makes the test
    exhaustive over the whole run instead of a sample.
    """
    seen: Dict[Tuple[Any, ...], Definition] = {}
    collisions: List[Dict[str, Any]] = []
    for definition in kept:
        found = seen.get(definition.signature)
        if found is None:
            seen[definition.signature] = definition
            continue
        if found.operator == definition.operator:
            continue
        collisions.append({
            "inventory_id": ion.inventory_id,
            "family_id": ion.family_id,
            "canonical_smiles": ion.canonical_smiles,
            "first_operator": found.operator, "first_mode": found.mode,
            "second_operator": definition.operator,
            "second_mode": definition.mode,
            "signature": definition.signature_text,
            "first_description": found.description,
            "second_description": definition.description})
    return collisions


def context_from(mol: Chem.Mol, core_atoms: Sequence[int],
                 inventory_id: int, family_id: str, smiles: str,
                 elements: Sequence[int]) -> Context:
    """A context over any molecule whose core atoms are known.

    Sequential editing needs this: after the first round the molecule being
    edited is a product, not a curated ion, and its core is known from the
    label certification already read rather than from a fresh search. The
    starting-ion path and this one share everything below the core, so an edit
    in round twelve is the same kind of edit as one in round one.
    """
    from .fragments import components_of, fixed_atoms
    core = frozenset(int(i) for i in core_atoms)
    fixed, _absorbed = fixed_atoms(mol, sorted(core))
    pieces: List[Component] = []
    for found in components_of(mol, sorted(core)):
        pieces.append(Component(frozenset(found.atom_indices),
                                found.attachment_roots[0],
                                found.core_neighbours[0],
                                found.attachment_orders[0], found.key))
    return Context(inventory_id=inventory_id, family_id=family_id,
                   canonical_smiles=smiles, mol=mol, core=core,
                   fixed=frozenset(fixed), components=tuple(pieces),
                   elements=tuple(elements))
