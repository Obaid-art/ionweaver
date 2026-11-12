"""Molecule depictions in the manuscript's own visual language.

The fills are the figures' fills, and they code chemistry roles only:

    blue      the protected core, carried through the edit untouched
    lavender  the neutral body outside the core, which an edit may change
    amber     what this edit added or changed
    rose      what it removed (drawn on the parent)

Atom labels stay black, so colour never means two things at once.  Every
molecule is drawn at one fixed bond length, which is what makes a grid of
structures read as a set rather than as differently zoomed pictures.  Output is
SVG on a transparent background.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

__all__ = ["molecule_svg", "ROLES", "CORE_FILL", "BODY_FILL", "NEW_FILL", "GONE_FILL"]

RGB = Tuple[float, float, float]
CORE_FILL: RGB = (0.62, 0.76, 0.92)   # blue
BODY_FILL: RGB = (0.85, 0.79, 0.93)   # lavender
NEW_FILL: RGB = (0.98, 0.78, 0.52)    # amber
GONE_FILL: RGB = (0.95, 0.72, 0.72)   # rose


def _hex(rgb: RGB) -> str:
    return "#" + "".join(f"{max(0, min(255, round(c * 255))):02X}" for c in rgb)


ROLES = [
    {"key": "core", "color": _hex(CORE_FILL), "label": "protected core"},
    {"key": "body", "color": _hex(BODY_FILL), "label": "body, outside the core"},
    {"key": "new", "color": _hex(NEW_FILL), "label": "added or changed by the edit"},
    {"key": "gone", "color": _hex(GONE_FILL), "label": "removed by the edit"},
]

BOND_PX = 30
_MIN_CANVAS = (200, 150)


def _median_bond(mol: Chem.Mol) -> float:
    conf = mol.GetConformer()
    lengths = []
    for b in mol.GetBonds():
        p = conf.GetAtomPosition(b.GetBeginAtomIdx())
        q = conf.GetAtomPosition(b.GetEndAtomIdx())
        lengths.append(((p.x - q.x) ** 2 + (p.y - q.y) ** 2) ** 0.5)
    return float(np.median(lengths)) if lengths else 1.5


def _canvas_for(mol: Chem.Mol, max_width: int, max_height: int) -> Tuple[int, int]:
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    xs = [conf.GetAtomPosition(i).x for i in range(n)]
    ys = [conf.GetAtomPosition(i).y for i in range(n)]
    scale = BOND_PX / _median_bond(mol)
    w = int((max(xs) - min(xs)) * scale * 1.16 + BOND_PX * 2.0)
    h = int((max(ys) - min(ys)) * scale * 1.16 + BOND_PX * 2.0)
    return (max(_MIN_CANVAS[0], min(int(max_width), w)),
            max(_MIN_CANVAS[1], min(int(max_height), h)))


def _prepared(mol: Chem.Mol) -> Chem.Mol:
    m = Chem.Mol(mol)
    for atom in m.GetAtoms():
        atom.SetAtomMapNum(0)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    return m


def molecule_svg(
    mol: Chem.Mol, *, core: Iterable[int] = (), added: Iterable[int] = (),
    removed: Iterable[int] = (), show_body: bool = True, q_labels: bool = False,
    max_width: int = 420, max_height: int = 300, notes: Optional[Dict[int, str]] = None,
    focus: Iterable[int] = (),
) -> str:
    """Draw one molecule with its chemistry roles washed in.

    ``q_labels`` annotates the core atoms Q0, Q1, ... in the order given;
    ``notes`` puts any text beside any atom (the studio numbers the body atoms);
    ``focus`` draws the named atoms in amber, on top of their role.
    """
    if mol is None or mol.GetNumAtoms() == 0:
        return empty_svg(min(max_width, 200), min(max_height, 150))
    m = _prepared(mol)
    n = m.GetNumAtoms()
    core_l = [int(i) for i in core if 0 <= int(i) < n]
    added_s = {int(i) for i in added if 0 <= int(i) < n} | {int(i) for i in focus if 0 <= int(i) < n}
    removed_s = {int(i) for i in removed if 0 <= int(i) < n}
    core_s = set(core_l) - {int(i) for i in focus}
    if q_labels:
        for q, idx in enumerate(core_l):
            m.GetAtomWithIdx(idx).SetProp("atomNote", f"Q{q}")
    for idx, text in (notes or {}).items():
        if 0 <= int(idx) < n:
            m.GetAtomWithIdx(int(idx)).SetProp("atomNote", str(text))

    fills: Dict[int, RGB] = {}
    for atom in m.GetAtoms():
        i = int(atom.GetIdx())
        if i in core_s:
            fills[i] = CORE_FILL
        elif i in added_s:
            fills[i] = NEW_FILL
        elif i in removed_s:
            fills[i] = GONE_FILL
        elif show_body:
            fills[i] = BODY_FILL
    bond_fills: Dict[int, RGB] = {}
    for bond in m.GetBonds():
        a = fills.get(int(bond.GetBeginAtomIdx()))
        b = fills.get(int(bond.GetEndAtomIdx()))
        if a is not None and a == b:
            bond_fills[int(bond.GetIdx())] = a

    w, h = _canvas_for(m, max_width, max_height)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    o = drawer.drawOptions()
    o.fixedBondLength = BOND_PX
    o.annotationFontScale = 0.85
    o.useBWAtomPalette()
    o.bondLineWidth = 2
    o.minFontSize = 10
    o.maxFontSize = 14
    o.padding = 0.06
    o.clearBackground = False
    o.highlightRadius = 0.36
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, m, highlightAtoms=list(fills), highlightAtomColors=fills,
        highlightBonds=list(bond_fills), highlightBondColors=bond_fills)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def empty_svg(width: int = 200, height: int = 150) -> str:
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}"><text x="50%" y="50%" text-anchor="middle" '
            f'font-family="system-ui,sans-serif" font-size="13" fill="#7C8896">'
            f'structure unavailable</text></svg>')
