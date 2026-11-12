"""Bring your own family: a sheet of SMILES and family labels.

The sheet goes through exactly the stages the published inventory went through,
with the engine's own functions:

1. the table is read (.csv, .tsv, .txt or .xlsx); the SMILES column and the
   family column are found by name, or by content when the headers are unusual;
   an ion-type or charge column is used when present and inferred from the
   structure when not;
2. ``curation.curate`` keeps one connected, closed-shell, metal-free monovalent
   ion per row and drops repeats within a family, recording why each row was
   kept or not;
3. ``registry.build`` learns the atom and bond environments the sheet displays;
4. ``discovery.discover`` finds each family's cores with the layered
   maximum-common-substructure search of the study (a family needs at least
   two retained ions);
5. ``fragments.starting_ions`` and ``fragments.build`` give the editable set
   and the sheet's own component library.

Edits on the sheet's ions are then decided against the published workspace
merged with the sheet's: the published library plus the sheet's components, and
the published observed environments plus the sheet's.
"""
from __future__ import annotations

import io
import re
import time
import uuid
from typing import Any, Dict, List, Optional, Sequence

import pandas as pd

from iw import curation, discovery, fragments, registry
from iw.common import CHEMISTRY_ERRORS, net_charge, parse

from .generate import Workspace, merged_workspace
from .library import Published

__all__ = ["read_table", "Sheet", "load_sheet", "STATUS_TEXT", "TEMPLATE_CSV", "NOVEL_EXAMPLE_CSV"]

TEMPLATE_CSV = ("SMILES,family\n"
                "CCCCn1cc[n+](C)c1,imidazolium\nCCCCCCn1cc[n+](C)c1,imidazolium\nCCn1cc[n+](C)c1,imidazolium\n"
                "CCCC[n+]1ccccc1,pyridinium\nCCCCCC[n+]1ccccc1,pyridinium\n"
                "CC(=O)[O-],acetates\nCCC(=O)[O-],acetates\nCCCC(=O)[O-],acetates\n")

#: a second example: families absent from the published inventory, plus the two cases a
#: user should expect to see reported rather than edited (a lone ion; shared part too small)
NOVEL_EXAMPLE_CSV = "\n".join([
    "SMILES,family",
    "Cc1cc(C)[o+]c(C)c1,pyrylium",
    "CCc1cc(C)[o+]c(C)c1,pyrylium",
    "CCCCc1cc(C)[o+]c(C)c1,pyrylium",
    "Cc1cc(-c2ccccc2)[o+]c(C)c1,pyrylium",
    "COc1ccc(-c2cc(C)[o+]c(C)c2)cc1,pyrylium",
    "C[S+](C)(C)=O,sulfoxonium",
    "CC[S+](C)(C)=O,sulfoxonium",
    "CCCC[S+](C)(C)=O,sulfoxonium",
    "CCCCCC[S+](C)(C)=O,sulfoxonium",
    "CN(C)C(=S)[S-],dithiocarbamate",
    "CCN(CC)C(=S)[S-],dithiocarbamate",
    "S=C([S-])N1CCCCC1,dithiocarbamate",
    "CCCCN(CCCC)C(=S)[S-],dithiocarbamate",
    "S=C([S-])N1CCOCC1,dithiocarbamate",
    "c1ccc([I+]c2ccccc2)cc1,iodonium",
    "Cc1ccc([I+]c2ccc(C)cc2)cc1,iodonium",
    "Cc1ccc([I+]c2ccccc2)cc1,iodonium",
    "O=C1[N-]S(=O)(=O)c2ccccc12,saccharinate",
    "CCCC[S-],thiolate",
    "c1ccc([S-])cc1,thiolate",
    "CC(C)(C)[S-],thiolate",
    "CCCCCCCC[S-],thiolate",
]) + "\n"

STATUS_TEXT = {
    "retained": "Kept",
    "excluded_out_of_scope": "Family outside the method's scope",
    "excluded_unreadable_structure": "SMILES could not be read",
    "excluded_disconnected": "More than one disconnected part",
    "excluded_radical": "Unpaired electron",
    "excluded_metal_atom": "Contains a metal atom",
    "excluded_unknown_ion_type": "Ion type is neither cation nor anion",
    "excluded_declared_charge_not_required": "Declared charge is not the ±1 its type requires",
    "excluded_computed_charge_not_required": "Net charge from the structure is not the ±1 its type requires",
    "excluded_charge_not_on_a_single_atom": "Charge is spread over more than one atom",
    "excluded_duplicate_in_family": "Repeats an ion already in this family",
    "missing_family": "No family label",
}

_SMILES_NAMES = ("smiles", "canonical smiles", "structure", "smile", "ion smiles")
_FAMILY_NAMES = ("family", "family label", "class", "ion family", "group", "label")
_ION_NAMES = ("ion type", "iontype", "ion", "type")
_CHARGE_NAMES = ("charge", "net charge", "formal charge", "declared charge")
#: a sheet larger than this is refused rather than left to run for minutes
MAX_ROWS = 3000
MAX_BYTES = 6 * 1024 * 1024
_MISSING = ("", "nan", "none", "null", "na", "n/a")


def _parses(text: str) -> bool:
    try:
        parse(text)
        return True
    except CHEMISTRY_ERRORS:
        return False


def read_table(data: bytes | str, filename: str = "") -> pd.DataFrame:
    name = (filename or "").lower()
    if len(data) > MAX_BYTES:
        raise ValueError("That file is larger than 6 MB. Split it, or keep the two columns the method needs.")
    if isinstance(data, bytes) and name.endswith((".xlsx", ".xlsm")):
        return pd.read_excel(io.BytesIO(data), dtype=str, keep_default_na=False)
    if isinstance(data, bytes) and name.endswith(".xls"):
        raise ValueError("Old .xls workbooks are not read; save the sheet as .xlsx or .csv.")
    text = data.decode("utf-8-sig", errors="replace") if isinstance(data, bytes) else str(data)
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        raise ValueError("That file is empty.")
    if len(lines) > MAX_ROWS + 1:
        raise ValueError(f"That file has more than {MAX_ROWS:,} rows; the shared-core search would take too long. "
                         "Split it by family.")
    head = lines[0]
    sep = next((c for c in ("\t", ",", ";") if c in head), None)
    # a file without a header line: its first field already reads as a molecule
    first = head.split(sep)[0].strip().strip('"') if sep else head.split()[0]
    headerless = _parses(first)
    if sep:
        if headerless:
            width = max(len(l.split(sep)) for l in lines[:50])
            names = ["SMILES", "family"] + [f"column {i}" for i in range(3, width + 1)]
            return pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, dtype=str, keep_default_na=False,
                               header=None, names=names[:width], index_col=False)
        return pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, dtype=str, keep_default_na=False, index_col=False)
    # whitespace separated: SMILES, then the rest of the line is the family label
    rows = [l.split(None, 1) for l in (lines if headerless else lines[1:])]
    cols = ["SMILES", "family"] if headerless else [c.strip() for c in head.split(None, 1)] + ["family"][len(head.split(None, 1)) - 1:]
    cols = (cols + ["family"])[:2]
    return pd.DataFrame([(r + [""])[:2] for r in rows], columns=cols)


def _column(df: pd.DataFrame, names: Sequence[str], *, exact: bool = False) -> Optional[Any]:
    """The column whose name is one of ``names``: exactly, or as whole words."""
    lookup = {re.sub(r"[\s_\-]+", " ", str(c).strip().lower()): c for c in df.columns}
    for n in names:
        if n in lookup:
            return lookup[n]
    if exact:
        return None
    for key, col in lookup.items():
        if any(re.search(r"(^|\s)%s(\s|$)" % re.escape(n), key) for n in names):
            return col
    return None


def _guess_smiles_column(df: pd.DataFrame) -> Optional[Any]:
    best, rate_best = None, 0.0
    for col in df.columns:
        values = [str(v) for v in df[col].head(30).tolist() if str(v).strip()]
        if not values:
            continue
        ok = 0
        for v in values:
            try:
                parse(v)
                ok += 1
            except CHEMISTRY_ERRORS:
                pass
        rate = ok / len(values)
        if rate > rate_best:
            best, rate_best = col, rate
    return best if rate_best >= 0.6 else None


def entries_from(df: pd.DataFrame) -> List[Dict[str, Any]]:
    """Rows of the sheet in the engine's entry format."""
    smiles_col = _column(df, _SMILES_NAMES) or _guess_smiles_column(df)
    if smiles_col is None:
        raise ValueError("No SMILES column was found. Name one column 'SMILES'.")
    family_col = _column(df, _FAMILY_NAMES)
    if family_col is None or family_col == smiles_col:
        others = [c for c in df.columns if c != smiles_col]
        family_col = others[0] if len(others) == 1 else None
    ion_col = _column(df, _ION_NAMES, exact=True)
    charge_col = _column(df, _CHARGE_NAMES, exact=True)
    if ion_col in (smiles_col, family_col):
        ion_col = None
    if charge_col in (smiles_col, family_col):
        charge_col = None

    if len(df) > MAX_ROWS:
        raise ValueError(f"That sheet has {len(df):,} rows; at most {MAX_ROWS:,} are read at once. Split it by family.")
    entries: List[Dict[str, Any]] = []
    smiles_v = df[smiles_col].astype(str).tolist()
    family_v = df[family_col].astype(str).tolist() if family_col is not None else [""] * len(df)
    ion_v = df[ion_col].astype(str).tolist() if ion_col is not None else [""] * len(df)
    charge_v = df[charge_col].astype(str).tolist() if charge_col is not None else [""] * len(df)
    clean = lambda x: "" if str(x).strip().lower() in _MISSING else str(x).strip()
    for k in range(len(df)):
        i = k + 2
        smi, fam, ion_type, charge = clean(smiles_v[k]), clean(family_v[k]), clean(ion_v[k]).lower(), clean(charge_v[k])
        if not smi:
            continue
        if ion_type not in ("cation", "anion") or not charge:
            try:
                net = net_charge(parse(smi))
            except CHEMISTRY_ERRORS:
                net = None
            if ion_type not in ("cation", "anion"):
                ion_type = "cation" if (net or 0) > 0 else ("anion" if (net or 0) < 0 else "unknown")
            if not charge and net is not None:
                charge = str(net)
        entries.append({
            "raw_index": len(entries) + 1, "source_workbook": "(uploaded)", "source_sheet": "",
            "source_row": i, "ion_type_raw": ion_type, "abbreviation_raw": "", "chemical_name_raw": "",
            "smiles_raw": smi, "declared_charge_raw": charge, "family_raw": fam,
        })
    if not entries:
        raise ValueError("No rows with a SMILES string were found.")
    return entries


class Sheet:
    """A user's family sheet after curation, core discovery and library building."""

    def __init__(self, P: Published, data: bytes | str, filename: str = "") -> None:
        t0 = time.perf_counter()
        self.sheet_id = uuid.uuid4().hex[:12]
        self.filename = filename or "pasted.csv"
        df = read_table(data, filename)
        entries = entries_from(df)
        # rows without a family label cannot be curated: say so instead of guessing
        labelled = [e for e in entries if e["family_raw"]]
        self.unlabelled = [e for e in entries if not e["family_raw"]]
        self.audit, self.inventory = curation.curate(labelled) if labelled else ([], [])
        for e in self.unlabelled:
            self.audit.append({**e, "family_id": "", "status": "missing_family", "structure_check": "",
                               "canonical_smiles": "", "heavy_atoms": ""})
        self.audit.sort(key=lambda r: r["source_row"])
        if not self.inventory:
            raise ValueError("No row survived curation. " + self._why_nothing())
        reg = registry.build(self.inventory)
        self.observed, self.observed_bonds = frozenset(reg["observed"]), frozenset(reg["observed_bonds"])
        self.elements = sorted({e[0] for e in reg["creatable"]})
        self.cores = discovery.discover(self.inventory)
        self.starting = fragments.starting_ions(self.cores, self.inventory)
        self.by_id = {r["inventory_id"]: r for r in self.starting}
        built = fragments.build(self.starting) if self.starting else {"library": [], "fragments": [], "instances": []}
        self.library, self.fragments = built["library"], built["fragments"]
        self.smiles_set = {r["canonical_smiles"] for r in self.inventory}
        self.workspace: Workspace = merged_workspace(
            P, f"sheet:{self.sheet_id}", self.library, self.fragments, self.observed, self.observed_bonds,
            self.elements, self.smiles_set)
        self.published_patterns = {m.pattern: m for m in P.core_models}
        published_keys = {e["component_key"] for e in P.library}
        self.components_new_to_published = sum(1 for e in self.library if e["component_key"] not in published_keys)
        self.elapsed = round(time.perf_counter() - t0, 2)

    def _why_nothing(self) -> str:
        from collections import Counter
        c = Counter(r["status"] for r in self.audit)
        return "; ".join(f"{STATUS_TEXT.get(k, k)}: {n}" for k, n in c.most_common())

    def record(self, inventory_id: int) -> Dict[str, Any]:
        try:
            return self.by_id[int(inventory_id)]
        except KeyError as exc:
            raise KeyError(f"no editable ion with id {inventory_id} in this sheet") from exc

    def ions_of(self, family_id: str) -> List[Dict[str, Any]]:
        return [r for r in self.starting if r["family_id"] == family_id]

    def as_dict(self) -> Dict[str, Any]:
        from collections import Counter
        members = Counter(r["family_id"] for r in self.inventory)
        labels = {r["family_id"]: r["family_label"] for r in self.inventory}
        families: Dict[str, Dict[str, Any]] = {}
        for fid in sorted(members):
            families[fid] = {"family_id": fid, "family_label": labels[fid],
                             "ion_type": "cation" if fid.startswith("cation") else "anion",
                             "members": members[fid], "editable": 0, "cores": [], "uncovered": []}
        for m in self.cores["models"]:
            pub = self.published_patterns.get(m["pattern"])
            families[m["family_id"]]["cores"].append({
                "grammar_id": m["grammar_id"], "pattern": m["pattern"], "atoms": m["atoms"], "bonds": m["bonds"],
                "support": m["support"], "support_fraction": m["support_fraction_of_family"],
                "copied_from": m["copied_from_inventory_id"],
                "published_equivalent": ({"family_label": pub.family_label, "grammar_id": pub.grammar_id}
                                         if pub else None),
                "members": [r["inventory_id"] for r in self.starting
                            if r["family_id"] == m["family_id"] and r["grammar_id"] == m["grammar_id"]]})
        for r in self.starting:
            families[r["family_id"]]["editable"] += 1
        for u in self.cores["uncovered"]:
            families[u["family_id"]]["uncovered"].append({"inventory_id": u["inventory_id"],
                                                          "smiles": u["canonical_smiles"], "reason": u["reason"]})
        status = Counter(r["status"] for r in self.audit)
        id_by_raw = {x["raw_index"]: x["inventory_id"] for x in self.inventory}
        return {
            "sheet_id": self.sheet_id, "filename": self.filename, "elapsed_s": self.elapsed,
            "rows": len(self.audit), "retained": len(self.inventory), "editable": len(self.starting),
            "families": list(families.values()),
            "status_counts": {k: {"n": n, "text": STATUS_TEXT.get(k, k)} for k, n in status.items()},
            "audit": [{"row": r["source_row"], "smiles": r["smiles_raw"], "family": r["family_raw"],
                       "ion_type": r["ion_type_raw"], "status": r["status"],
                       "text": STATUS_TEXT.get(r["status"], r["status"]),
                       "canonical_smiles": r.get("canonical_smiles", ""),
                       "inventory_id": id_by_raw.get(r.get("raw_index")) if r["status"] == "retained" else None}
                      for r in self.audit],
            "ions": [{k: r[k] for k in ("inventory_id", "family_id", "family_label", "canonical_smiles", "grammar_id",
                                        "core_pattern", "heavy_atoms", "body_heavy_atoms")}
                     | {"core_atoms": [int(i) for i in str(r["core_atoms"]).split()]} for r in self.starting],
            "library_components": len(self.library), "library_new_to_published": self.components_new_to_published,
            "atom_environments": len(self.observed), "bond_environments": len(self.observed_bonds),
        }


_SHEETS: Dict[str, Sheet] = {}


def load_sheet(P: Published, data: bytes | str, filename: str = "") -> Sheet:
    sheet = Sheet(P, data, filename)
    if len(_SHEETS) > 40:
        oldest = next(iter(_SHEETS))
        _SHEETS.pop(oldest, None)
    _SHEETS[sheet.sheet_id] = sheet
    return sheet


def get_sheet(sheet_id: str) -> Optional[Sheet]:
    return _SHEETS.get(str(sheet_id))
