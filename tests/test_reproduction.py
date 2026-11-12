"""The server reproduces the submitted study.

Run from the repository root:  pytest -q
"""
from __future__ import annotations

import collections
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "server"))

from ionweaver_web import engine_fingerprint  # noqa: E402
from ionweaver_web.generate import (OPERATORS, capacity, chain, definitions_for, generate,  # noqa: E402
                                    published_workspace)
from ionweaver_web.identify import identify  # noqa: E402
from ionweaver_web.library import SEED, load_published  # noqa: E402
from ionweaver_web.sheet import TEMPLATE_CSV, load_sheet  # noqa: E402


@pytest.fixture(scope="session")
def P():
    return load_published()


@pytest.fixture(scope="session")
def WS(P):
    return published_workspace(P)


def test_engine_is_the_notebooks():
    fp = engine_fingerprint()
    assert fp["unmodified"], fp
    assert fp["pins"]["RDKIT_VERSION"] == "2026.03.6"


def test_library_reproduces_the_archive(P):
    assert P.checks["all"], P.checks
    c = P.counts()
    assert (c["ions"], c["families"], c["cores"], c["editable_ions"], c["library_components"]) == (1241, 38, 44, 1222, 255)
    assert (c["atom_environments"], c["bond_environments"], c["terminal_fragments"]) == (57, 71, 247)


def test_capacities_match_the_published_definition_counts(P, WS):
    summary = {(int(r["inventory_id"]), r["operator"]): int(r["definitions_available"]) for r in P.per_ion_operator}
    for iid in (1, 2, 300, 455, 458, 700, 900, 1088, 1089, 1200):
        _ion, kept, _tally = definitions_for(P.record(iid), WS)
        caps = capacity(kept)
        for op in OPERATORS:
            if (iid, op) in summary:
                assert caps[op] == summary[(iid, op)], (iid, op)


def test_products_match_the_published_examination_order(P, WS):
    rows = list(csv.DictReader(open(ROOT / "tests" / "fixtures" / "examination_records_sample.csv",
                                    newline="", encoding="utf-8")))
    by = collections.defaultdict(list)
    for r in rows:
        by[int(r["inventory_id"])].append(r)
    for iid, pub in by.items():
        rec = P.record(iid)
        res = generate(rec, WS, want=10 ** 6, seed=SEED, svg=False, max_examinations=10 ** 6)
        mine = {(p["operator"], p["order_index"]): p["product_smiles"] for p in res["products"]}
        delivered = set(mine.values())
        checked = 0
        for r in pub:
            key = (r["operator"], int(r["order_index"]))
            if r["certified"] == "1" and r["identical_to_parent"] == "0" and r["repeat_for_this_ion"] == "0":
                # the same product at the same position, or the same product already
                # delivered under another operator (the study ran operators one after
                # another, this server interleaves them, so which one "owns" a product
                # two operators both reach can differ; the product set cannot)
                assert mine.get(key) == r["product_smiles"] or r["product_smiles"] in delivered, (iid, key)
                checked += 1
        assert checked > 0
        assert res["summary"]["refused_by_certification"] == 0


def test_every_delivered_product_is_distinct_and_certified(P, WS):
    res = generate(P.record(458), WS, want=200, svg=False)
    smiles = [p["product_smiles"] for p in res["products"]]
    assert len(smiles) == len(set(smiles)) == 200
    assert all(p["certified"] and all(p["checks"].values()) for p in res["products"])
    assert all(p["judge_agrees"] for p in res["products"])
    assert P.record(458)["canonical_smiles"] not in smiles


def test_identify_recognises_reference_ions(P):
    v = identify("CCCCn1cc[n+](C)c1", P)
    assert v["ok"] and v["best"]["family_label"] == "imidazolium" and v["published_id"] == 42
    assert not identify("CCOCC", P)["ok"]
    assert not identify("[Na+]", P)["ok"]


def test_example_sheet_discovers_cores_and_edits(P):
    sheet = load_sheet(P, TEMPLATE_CSV, "example.csv")
    d = sheet.as_dict()
    assert d["retained"] == 8 and d["editable"] == 8 and len(d["families"]) == 3
    res = generate(sheet.record(1), sheet.workspace, want=10, svg=False)
    assert res["summary"]["delivered"] == 10


def test_studio_moves_agree_with_the_generator(P, WS):
    from ionweaver_web.studio import catalogue, state_from_record, try_move, move_key, _definitions, _context
    rec = P.record(42)
    st = state_from_record(rec, WS)
    cat = catalogue(st, WS, svg=False)
    assert cat["total"] == sum(cat["capacity"].values())
    # every kept definition has a key the catalogue lists, and every listed move applies
    ion = _context(st, WS)
    kept = _definitions(ion, WS)
    keys = {move_key(d) for d in kept}
    assert len(keys) == len(kept)                      # no two definitions share a user-facing name
    ok = try_move(st, WS, {"operator": "Graft", "atom": 0, "structure": "*CCO"}, svg=False)
    assert ok["ok"] and ok["record"]["certified"] and len(ok["state"]["core_atoms"]) == 5
    # the product the studio built is the product the generator certifies: re-run its catalogue
    again = catalogue(ok["state"], WS, svg=False)
    assert again["total"] > 0
    # a refused move names the rule
    core_atom = st["core_atoms"][0]
    no = try_move(st, WS, {"operator": "Substitute", "atom": core_atom, "element": "S"}, svg=False)
    assert not no["ok"] and "core" in no["condition"]
    no2 = try_move(st, WS, {"operator": "Graft", "atom": 0, "structure": "*C"}, svg=False)
    assert not no2["ok"] and no2["condition"] == "one atom"


NOVEL_SHEET = """SMILES,family
Cc1cc(C)[o+]c(C)c1,pyrylium
CCc1cc(C)[o+]c(C)c1,pyrylium
CCCCc1cc(C)[o+]c(C)c1,pyrylium
Cc1cc(-c2ccccc2)[o+]c(C)c1,pyrylium
C[S+](C)(C)=O,sulfoxonium
CC[S+](C)(C)=O,sulfoxonium
CCCC[S+](C)(C)=O,sulfoxonium
CN(C)C(=S)[S-],dithiocarbamate
CCN(CC)C(=S)[S-],dithiocarbamate
S=C([S-])N1CCCCC1,dithiocarbamate
c1ccc([I+]c2ccccc2)cc1,iodonium
Cc1ccc([I+]c2ccc(C)cc2)cc1,iodonium
O=C1[N-]S(=O)(=O)c2ccccc12,saccharinate
"""


def test_families_nobody_published_get_their_own_core(P):
    from rdkit import Chem
    sheet = load_sheet(P, NOVEL_SHEET, "novel.csv")
    d = sheet.as_dict()
    found = {f["family_label"]: f for f in d["families"]}
    for label, charged in (("pyrylium", "[#8+]"), ("sulfoxonium", "[#16+]"), ("dithiocarbamate", "[#16-]"), ("iodonium", "[#53+]")):
        fam = found[label]
        assert fam["editable"] == fam["members"] and len(fam["cores"]) == 1, label
        assert charged in fam["cores"][0]["pattern"] and fam["cores"][0]["published_equivalent"] is None, label
    assert found["saccharinate"]["editable"] == 0 and found["saccharinate"]["uncovered"]
    # every product of a new family keeps the new core and the charge
    for rec in [next(r for r in sheet.starting if r["family_label"] == lab) for lab in ("pyrylium", "dithiocarbamate", "iodonium")]:
        res = generate(rec, sheet.workspace, want=25, svg=False)
        q = Chem.MolFromSmarts(rec["core_pattern"])
        assert res["summary"]["delivered"] == 25 and res["summary"]["refused_by_certification"] == 0
        for p in res["products"]:
            m = Chem.MolFromSmiles(p["product_smiles"])
            assert m.HasSubstructMatch(q) and all(p["checks"].values())
    # a novel ion on its own is told to bring its family
    assert identify("Cc1cc(C)[o+]c(C)c1", P)["status"] == "no_core"


def test_body_rules_shape_what_the_edits_build(P, WS):
    from ionweaver_web.generate import batch
    rules = {"body_elements": ["C", "N", "O"], "body_multiple_bonds": [None, 0], "body_rings": [None, 0], "hbd": [None, 1]}
    res = generate(P.record(42), WS, want=30, constraints=rules, svg=False)
    assert res["summary"]["delivered"] == 30
    for p in res["products"]:
        d = p["descriptors"]
        assert set(d["body_elements"]) <= {"C", "N", "O"} and d["body_multiple_bonds"] == 0 and d["body_rings"] == 0 and d["hbd"] <= 1
    # the same rules across families, in one batch; a family whose parents cannot comply says so
    b = batch(P, WS, family_ids=["cation__pyridinium", "anion__bf4_derivatives"], per_family=10, ions_per_family=3,
              constraints={"body_elements": ["C"], "body_multiple_bonds": [None, 0]}, svg=False)
    rows = {f["family_label"]: f for f in b["families"]}
    assert rows["pyridinium"]["delivered"] == 10
    assert rows["BF4 derivatives"]["parents_breaking_rules"] == 3      # every parent's body is fluorine
    assert all(set(p["descriptors"]["body_elements"]) <= {"C"} for p in b["products"])


def test_chain_moves_and_keeps_the_core(P, WS):
    out = chain(P.record(458), WS, rule="exploration", depth=4, svg=False)
    assert out["summary"]["rounds_moved"] == 4
    assert all(all(s["selected"]["checks"].values()) for s in out["steps"] if s.get("moved"))
