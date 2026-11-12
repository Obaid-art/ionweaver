"""Quick RDKit descriptors for filtering a set of generated ions.

Every value is computed by RDKit's own descriptor functions on the product as
generated.  They are properties of one ion, not of an ionic liquid, and they are
offered as filters, never as a claim about any property of a salt.
"""
from __future__ import annotations

from typing import Any, Dict, Iterable, List

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

__all__ = ["DESCRIPTORS", "BODY_RULES", "descriptors_of", "body_profile", "ELEMENT_ORDER"]

#: key, label, unit, how it is computed, decimals
DESCRIPTORS: List[Dict[str, Any]] = [
    {"key": "mw", "label": "Molecular weight", "unit": "g/mol", "decimals": 1,
     "how": "RDKit Descriptors.MolWt (average atomic masses)"},
    {"key": "heavy_atoms", "label": "Atoms (H not counted)", "unit": "", "decimals": 0,
     "how": "atoms other than hydrogen"},
    {"key": "hba", "label": "H-bond acceptors", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumHBA"},
    {"key": "hbd", "label": "H-bond donors", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumHBD"},
    {"key": "logp", "label": "Crippen log P", "unit": "", "decimals": 2,
     "how": "Crippen.MolLogP; an atom-contribution estimate on a charged ion, use with care"},
    {"key": "tpsa", "label": "Polar surface area", "unit": "Å²", "decimals": 1,
     "how": "rdMolDescriptors.CalcTPSA"},
    {"key": "rotatable", "label": "Rotatable bonds", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumRotatableBonds"},
    {"key": "rings", "label": "Rings", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumRings"},
    {"key": "aromatic_rings", "label": "Aromatic rings", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumAromaticRings"},
    {"key": "fsp3", "label": "Fraction sp³ carbon", "unit": "", "decimals": 2,
     "how": "rdMolDescriptors.CalcFractionCSP3"},
    {"key": "heteroatoms", "label": "Heteroatoms", "unit": "", "decimals": 0,
     "how": "rdMolDescriptors.CalcNumHeteroatoms"},
    {"key": "halogens", "label": "Halogen atoms", "unit": "", "decimals": 0,
     "how": "count of F, Cl, Br, I"},
    # -- the body: everything outside the protected core (and outside any ring it holds) --
    {"key": "body_heavy_atoms", "label": "Atoms outside the core (H not counted)", "unit": "", "decimals": 0, "body": True,
     "how": "atoms not in the protected core nor in a ring the core holds"},
    {"key": "body_multiple_bonds", "label": "Double or triple bonds outside the core", "unit": "", "decimals": 0, "body": True,
     "how": "non-aromatic double and triple bonds with at least one end outside the core"},
    {"key": "body_rings", "label": "Rings outside the core", "unit": "", "decimals": 0, "body": True,
     "how": "rings not entirely inside the protected core"},
    {"key": "body_aromatic_rings", "label": "Aromatic rings outside the core", "unit": "", "decimals": 0, "body": True,
     "how": "aromatic rings not entirely inside the protected core"},
    {"key": "body_heteroatoms", "label": "Heteroatoms outside the core", "unit": "", "decimals": 0, "body": True,
     "how": "non-carbon atoms outside the core, hydrogen not counted"},
]

#: the rules a user may put on the body, besides numeric bounds on the rows above
BODY_RULES = [
    {"key": "body_elements", "label": "Elements allowed outside the core",
     "how": "every atom outside the core must be one of the allowed elements"},
]

ELEMENT_ORDER = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Si", "B", "Se"]

_HALOGENS = {9, 17, 35, 53}


def descriptors_of(mol: Chem.Mol) -> Dict[str, Any]:
    """Every descriptor in :data:`DESCRIPTORS`, plus the molecular formula."""
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(float(Descriptors.MolWt(mol)), 2),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "hba": int(rdMolDescriptors.CalcNumHBA(mol)),
        "hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 2),
        "rotatable": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "rings": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fsp3": round(float(rdMolDescriptors.CalcFractionCSP3(mol)), 3),
        "heteroatoms": int(rdMolDescriptors.CalcNumHeteroatoms(mol)),
        "halogens": sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _HALOGENS),
    }


def body_profile(mol: Chem.Mol, core: Iterable[int]) -> Dict[str, Any]:
    """What lies outside the protected core: elements, multiple bonds, rings, size.

    The core and any ring it holds in part are "fixed"; the body is the rest.
    """
    from iw.fragments import fixed_atoms
    core_l = [int(i) for i in core]
    fixed, _absorbed = fixed_atoms(mol, core_l)
    fixed = set(fixed)
    body = [a for a in mol.GetAtoms() if int(a.GetIdx()) not in fixed]
    multiple = 0
    for b in mol.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in fixed and j in fixed:
            continue
        if b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
            multiple += 1
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    outside = [r for r in rings if not r <= fixed]
    aromatic = [r for r in outside if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    return {
        "body_heavy_atoms": len(body),
        "body_elements": sorted({a.GetSymbol() for a in body}, key=lambda e: (ELEMENT_ORDER.index(e) if e in ELEMENT_ORDER else 99, e)),
        "body_multiple_bonds": multiple,
        "body_rings": len(outside),
        "body_aromatic_rings": len(aromatic),
        "body_heteroatoms": sum(1 for a in body if a.GetAtomicNum() not in (1, 6)),
    }
