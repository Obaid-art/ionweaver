
# -*- coding: utf-8 -*-
"""Stage 05b — realising a definition, and certifying what came out.

The definition was decided in `operators` from the parent graph alone. Nothing
here decides anything: the builder's only job is to produce the graph the
definition specified, in the state the ledger predicted, and then to say so.

Every atom is labelled before the edit and read after it. The core's atoms
carry a label applied before anything changes, so requirements 3 and 4 are
tested on the atoms the core first matched and no fresh substructure search is
ever run on a product. The label also gives the atom correspondence that the
mapped product SMILES carries, which is what lets a reader check requirement 4
without re-deriving the match.

Certification is an assertion and not a filter. It runs on every product as
a standing check that the construction rules held. Anything it ever refuses is
recorded with the requirement and the reason, and across a run it is expected to
refuse nothing; that expectation is the claim, not a pass rate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Set, Tuple

from rdkit import Chem

from .common import RunError, atom_environment
from .ledger import BOND_TYPE, Resulting
from .operators import Definition

PARENT = "iw_parent"
CORE = "iw_core"
#: An atom an edit introduced carries the structure it came from and its index
#: within it, so the perception check of C2 can put every ring and every
#: aromatic flag in a product beside the one its source carried.
ADDED = "iw_added"
SOURCE = "iw_source"

#: The five requirements of C3, in the order they are tested. The order matters
#: only for which one is named when a product is refused.
REQUIREMENTS: Tuple[str, ...] = (
    "1_one_connected_closed_shell_molecule_that_sanitizes",
    "2_the_formal_charge_its_family_requires",
    "3_every_core_atom_present_as_it_was",
    "4_the_induced_bond_set_on_those_atoms_unchanged",
    "5_all_formal_charge_inside_the_core",
)


@dataclass
class Product:
    """What a definition produced, and how it was produced."""
    mol: Optional[Chem.Mol] = None
    smiles: str = ""
    mapped_smiles: str = ""
    hydrogens_written_explicitly: Tuple[int, ...] = ()
    built: bool = False
    failure: str = ""
    reason: str = ""


def label(parent: Chem.Mol, core_atoms: Sequence[int]) -> Chem.Mol:
    """Copy the parent with every atom labelled, and the core's atoms marked."""
    work = Chem.Mol(parent)
    core = set(int(i) for i in core_atoms)
    for atom in work.GetAtoms():
        index = int(atom.GetIdx())
        atom.SetIntProp(PARENT, index)
        atom.SetIntProp(CORE, 1 if index in core else 0)
    return work


def _attach(work: Chem.RWMol, structure: str,
            anchors: Sequence[Tuple[int, str]]) -> List[int]:
    """Join one structure on, at the atoms and orders the definition names.

    The structure is written with a dummy standing at each attachment. A dummy
    says which of the structure's atoms is a root there; it is not added to the
    product. A library component carries one dummy. The only structure carrying
    two is the single atom a bridged insertion puts into a bond, where both
    dummies stand on that one atom, so pairing the dummies with the anchors
    cannot be ambiguous. Anything else with two would be, and stops the run.
    """
    piece = Chem.MolFromSmiles(structure)
    if piece is None:
        raise RunError("a definition names the structure %r, which does not "
                       "read as a molecule" % structure)
    dummies = [int(a.GetIdx()) for a in piece.GetAtoms()
               if a.GetAtomicNum() == 0]
    if len(dummies) != len(anchors):
        raise RunError("the structure %r carries %d attachment points and the "
                       "definition names %d anchors"
                       % (structure, len(dummies), len(anchors)))
    heavy = piece.GetNumAtoms() - len(dummies)
    if len(dummies) > 1 and heavy != 1:
        raise RunError("the structure %r attaches at %d points and holds %d "
                       "heavy atoms. Only a one-atom bridge attaches at more "
                       "than one point, where the pairing of attachments to "
                       "anchors cannot matter." % (structure, len(dummies),
                                                   heavy))

    mapping: Dict[int, int] = {}
    for atom in piece.GetAtoms():
        if int(atom.GetIdx()) in dummies:
            continue
        index = work.AddAtom(Chem.Atom(atom))
        mapping[int(atom.GetIdx())] = index
        placed_atom = work.GetAtomWithIdx(index)
        placed_atom.SetIntProp(ADDED, int(atom.GetIdx()))
        placed_atom.SetProp(SOURCE, structure)
    for bond in piece.GetBonds():
        a, b = int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx())
        if a in dummies or b in dummies:
            continue
        new = work.AddBond(mapping[a], mapping[b], bond.GetBondType()) - 1
        work.GetBondWithIdx(new).SetIsAromatic(bond.GetIsAromatic())

    for dummy, (anchor, order) in zip(dummies, anchors):
        root = int(piece.GetAtomWithIdx(dummy).GetBonds()[0]
                   .GetOtherAtomIdx(dummy))
        work.AddBond(mapping[root], int(anchor), BOND_TYPE[order])
    return sorted(mapping.values())


def realise(labelled: Chem.Mol, definition: Definition) -> Product:
    """Produce the graph the definition specified, in the predicted state.

    The hydrogen counts the ledger computed are the authority. They are left
    implicit where RDKit arrives at them on its own, which keeps the product's
    canonical form the ordinary one, and written explicitly where it does not.
    An atom whose count cannot be made to agree stops the run, because the
    condition that admitted the edit was then decided on a state the product is
    not in.
    """
    work = Chem.RWMol(labelled)

    for index, element in definition.element_changes:
        work.GetAtomWithIdx(int(index)).SetAtomicNum(int(element))
    for a, b, order in definition.order_changes:
        bond = work.GetBondBetweenAtoms(int(a), int(b))
        if bond is None:
            raise RunError("a definition changes the order of the bond %d-%d, "
                           "which the parent does not carry" % (a, b))
        bond.SetBondType(BOND_TYPE[order])
    for a, b in definition.remove_bonds:
        if work.GetBondBetweenAtoms(int(a), int(b)) is None:
            raise RunError("a definition removes the bond %d-%d, which the "
                           "parent does not carry" % (a, b))
        work.RemoveBond(int(a), int(b))
    for a, b, order in definition.add_bonds:
        work.AddBond(int(a), int(b), BOND_TYPE[order])

    for attachment in definition.attachments:
        _attach(work, attachment.structure, attachment.anchors)

    for index in sorted((int(i) for i in definition.remove_atoms),
                        reverse=True):
        work.RemoveAtom(index)

    mol = work.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except (ValueError, Chem.MolSanitizeException):
        mol = _write_hydrogens(work, definition, mol)
        if mol is None:
            return Product(failure="1_one_connected_closed_shell_molecule_"
                                   "that_sanitizes",
                           reason="the product does not sanitize even with the "
                                  "hydrogen counts the ledger computed")
        explicit = tuple(sorted(p.index for p in definition.results))
        return _finish(mol, definition, explicit)

    if _hydrogens_agree(mol, definition):
        return _finish(mol, definition, ())

    rewritten = _write_hydrogens(work, definition, mol)
    if rewritten is None:
        return Product(failure="1_one_connected_closed_shell_molecule_that_"
                               "sanitizes",
                       reason="the hydrogen counts the ledger computed cannot "
                              "be written onto the product")
    return _finish(rewritten, definition,
                   tuple(sorted(p.index for p in definition.results)))


def _where(mol: Chem.Mol) -> Dict[int, int]:
    """Parent atom index to the index it now sits at."""
    return {int(a.GetIntProp(PARENT)): int(a.GetIdx())
            for a in mol.GetAtoms() if a.HasProp(PARENT)}


def _hydrogens_agree(mol: Chem.Mol, definition: Definition) -> bool:
    """Did RDKit arrive at the counts the ledger computed on its own?"""
    where = _where(mol)
    for result in definition.results:
        if result.index < 0:
            continue                      # an atom the edit built
        index = where.get(result.index)
        if index is None:
            return False
        if int(mol.GetAtomWithIdx(index).GetTotalNumHs()) != result.hydrogens:
            return False
    return True


def _write_hydrogens(work: Chem.RWMol, definition: Definition,
                     mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Write the ledger's counts onto the touched atoms and sanitize again."""
    again = Chem.RWMol(work)
    where = _where(again)
    for result in definition.results:
        if result.index < 0:
            continue
        index = where.get(result.index)
        if index is None:
            continue
        atom = again.GetAtomWithIdx(index)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(int(result.hydrogens))
    out = again.GetMol()
    try:
        Chem.SanitizeMol(out)
    except (ValueError, Chem.MolSanitizeException):
        return None
    del mol
    return out


def _finish(mol: Chem.Mol, definition: Definition,
            explicit: Tuple[int, ...]) -> Product:
    """Canonicalise, and write the mapped form beside it."""
    smiles = Chem.MolToSmiles(mol)
    mapped = Chem.Mol(mol)
    for atom in mapped.GetAtoms():
        if atom.HasProp(PARENT):
            atom.SetAtomMapNum(int(atom.GetIntProp(PARENT)) + 1)
    return Product(mol=mol, smiles=smiles,
                   mapped_smiles=Chem.MolToSmiles(mapped),
                   hydrogens_written_explicitly=explicit, built=True)


# ---------------------------------------------------------------------------
# certification


@dataclass
class Certificate:
    """What certification found, requirement by requirement."""
    certified: bool = False
    failed: str = ""
    reason: str = ""
    core_atoms_found: int = 0
    core_bonds_found: int = 0
    core_bonds_expected: int = 0
    net_charge: int = 0
    charge_outside_the_core: int = 0
    fragments: int = 0
    radical_electrons: int = 0


def core_induced_bonds(mol: Chem.Mol,
                       labelled_core: FrozenSet[int]) -> Set[Tuple[int, int, str]]:
    """The complete induced bond set on the core's atoms, keyed by label.

    Both endpoints must be core atoms, which is what makes this the *induced*
    set: a bond from the core out to the body is not in it, and is the reason
    Rebond may change such a bond without touching the core.
    """
    out: Set[Tuple[int, int, str]] = set()
    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        if not (begin.HasProp(PARENT) and end.HasProp(PARENT)):
            continue
        a, b = int(begin.GetIntProp(PARENT)), int(end.GetIntProp(PARENT))
        if a not in labelled_core or b not in labelled_core:
            continue
        low, high = (a, b) if a < b else (b, a)
        out.add((low, high, str(bond.GetBondType())))
    return out


def certify(product: Chem.Mol, parent: Chem.Mol, required_charge: int,
            expected_bonds: Set[Tuple[int, int, str]],
            core_size: int) -> Certificate:
    """The five requirements of C3, as an assertion over a finished product.

    Requirements 3 and 4 read the label applied before the edit, so no fresh
    substructure search is run here. That is the whole point of labelling: a
    search on the product could find the core somewhere else, and then a
    product that moved the core would certify.
    """
    found = Certificate()
    found.fragments = len(Chem.GetMolFrags(product))
    found.radical_electrons = sum(int(a.GetNumRadicalElectrons())
                                  for a in product.GetAtoms())
    found.net_charge = sum(int(a.GetFormalCharge()) for a in product.GetAtoms())

    core_atoms = [a for a in product.GetAtoms()
                  if a.HasProp(CORE) and int(a.GetIntProp(CORE)) == 1]
    found.core_atoms_found = len(core_atoms)
    found.charge_outside_the_core = sum(
        int(a.GetFormalCharge()) for a in product.GetAtoms()
        if not (a.HasProp(CORE) and int(a.GetIntProp(CORE)) == 1))

    if found.fragments != 1 or found.radical_electrons:
        found.failed = REQUIREMENTS[0]
        found.reason = ("%d disconnected parts and %d unpaired electrons"
                        % (found.fragments, found.radical_electrons))
        return found

    if found.net_charge != required_charge:
        found.failed = REQUIREMENTS[1]
        found.reason = ("net charge %+d where the family requires %+d"
                        % (found.net_charge, required_charge))
        return found

    if found.core_atoms_found != core_size:
        found.failed = REQUIREMENTS[2]
        found.reason = ("%d of the core's %d atoms survive the edit"
                        % (found.core_atoms_found, core_size))
        return found
    labels = frozenset(int(a.GetIntProp(PARENT)) for a in core_atoms)
    for atom in core_atoms:
        # The label on a parent atom is that atom's own index, by construction
        # in `label`, so the parent is indexed directly and not searched.
        was = parent.GetAtomWithIdx(int(atom.GetIntProp(PARENT)))
        if (int(atom.GetAtomicNum()) != int(was.GetAtomicNum())
                or int(atom.GetFormalCharge()) != int(was.GetFormalCharge())
                or bool(atom.GetIsAromatic()) != bool(was.GetIsAromatic())):
            found.failed = REQUIREMENTS[2]
            found.reason = ("core atom %d is %s%+d%s in the product and "
                            "%s%+d%s in the parent"
                            % (int(atom.GetIntProp(PARENT)),
                               atom.GetSymbol(), atom.GetFormalCharge(),
                               " aromatic" if atom.GetIsAromatic() else "",
                               was.GetSymbol(), was.GetFormalCharge(),
                               " aromatic" if was.GetIsAromatic() else ""))
            return found

    actual = core_induced_bonds(product, labels)
    found.core_bonds_found = len(actual)
    found.core_bonds_expected = len(expected_bonds)
    if actual != expected_bonds:
        missing = sorted(expected_bonds - actual)
        extra = sorted(actual - expected_bonds)
        found.failed = REQUIREMENTS[3]
        found.reason = ("the induced bond set on the core differs: %d missing "
                        "%s, %d added %s"
                        % (len(missing), missing[:3], len(extra), extra[:3]))
        return found

    if found.charge_outside_the_core != 0:
        found.failed = REQUIREMENTS[4]
        found.reason = ("%+d of formal charge sits outside the core"
                        % found.charge_outside_the_core)
        return found

    found.certified = True
    return found
