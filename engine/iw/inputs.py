
# -*- coding: utf-8 -*-
"""Stage 01: reading the two property compilations.

Every column is resolved from a declared map rather than by position, and each
map carries an assertion against a header cell the sheet names. A map that does
not verify stops the run.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence

import openpyxl

from .common import RunError, sha256_file


@dataclass(frozen=True)
class SheetMap:
    """How one ion sheet is read. Positions are declared, not discovered."""
    workbook: str
    sheet: str
    ion_type: int
    smiles: int
    charge: int
    family: int
    abbreviation: Optional[int] = None
    chemical_name: Optional[int] = None
    #: header cells that must appear where stated, as the proof of the map
    assertions: Sequence[tuple] = ()


#: The two sheets the study starts from. Their layouts differ: the density
#: workbook names every column, the viscosity workbook names only two, so the
#: charge sits one column earlier there. Both are stated rather than inferred.
SHEET_MAPS: Sequence[SheetMap] = (
    SheetMap(workbook="density_P.xlsx", sheet="ion",
             ion_type=0, abbreviation=1, chemical_name=2, smiles=3,
             charge=5, family=6,
             assertions=((0, "Ion type"), (3, "SMILES"), (5, "Charge"),
                         (6, "Family"))),
    SheetMap(workbook="viscosity_P.xlsx", sheet="ions",
             ion_type=0, abbreviation=1, chemical_name=2, smiles=3,
             charge=4, family=5,
             assertions=((3, "SMILES"), (5, "Family"))),
)


def _norm(cell: Any) -> str:
    return "" if cell is None else str(cell).strip()


def read_sheet(directory: str, spec: SheetMap) -> List[Dict[str, Any]]:
    """Every row of one ion sheet, with its provenance."""
    path = os.path.join(directory, spec.workbook)
    if not os.path.exists(path):
        raise RunError("input workbook not found: %s" % path)

    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    if spec.sheet not in book.sheetnames:
        raise RunError("%s has no sheet %r (has %s)"
                       % (spec.workbook, spec.sheet, book.sheetnames))
    sheet = book[spec.sheet]
    rows = list(sheet.iter_rows(values_only=True))
    book.close()
    if not rows:
        raise RunError("%s!%s is empty" % (spec.workbook, spec.sheet))

    header = rows[0]
    for index, expected in spec.assertions:
        got = _norm(header[index] if index < len(header) else None)
        if got.lower() != expected.lower():
            raise RunError(
                "%s!%s column %d is %r, the map says %r. The map is wrong or "
                "the workbook changed; nothing is read positionally."
                % (spec.workbook, spec.sheet, index, got, expected))

    out: List[Dict[str, Any]] = []
    for number, row in enumerate(rows[1:], start=2):
        def at(idx: Optional[int]) -> str:
            if idx is None or idx >= len(row):
                return ""
            return _norm(row[idx])
        if not any(_norm(c) for c in row):
            continue
        out.append({
            "source_workbook": spec.workbook,
            "source_sheet": spec.sheet,
            "source_row": number,
            "ion_type_raw": at(spec.ion_type),
            "abbreviation_raw": at(spec.abbreviation),
            "chemical_name_raw": at(spec.chemical_name),
            "smiles_raw": at(spec.smiles),
            "declared_charge_raw": at(spec.charge),
            "family_raw": at(spec.family),
        })
    return out


def read_all(directory: str) -> List[Dict[str, Any]]:
    """Both sheets, in the order the maps are declared, with a raw index."""
    entries: List[Dict[str, Any]] = []
    for spec in SHEET_MAPS:
        entries.extend(read_sheet(directory, spec))
    for index, entry in enumerate(entries, start=1):
        entry["raw_index"] = index
    return entries


def provenance(directory: str) -> Dict[str, Any]:
    """What was read, and the hash of each file it was read from."""
    files = {}
    for spec in SHEET_MAPS:
        path = os.path.join(directory, spec.workbook)
        files[spec.workbook] = {
            "sha256": sha256_file(path),
            "bytes": os.path.getsize(path),
            "sheet_read": spec.sheet,
            "column_map": {
                "ion_type": spec.ion_type,
                "abbreviation": spec.abbreviation,
                "chemical_name": spec.chemical_name,
                "smiles": spec.smiles,
                "declared_charge": spec.charge,
                "family": spec.family,
            },
            "header_assertions": [{"column": i, "expected": e}
                                  for i, e in spec.assertions],
        }
    return {"source_files": files}
