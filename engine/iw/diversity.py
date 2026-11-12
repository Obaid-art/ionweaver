
# -*- coding: utf-8 -*-
"""Stage 10b -- what the generated set looks like, measured inside the run.

Validity, uniqueness, novelty and internal diversity are the four quantities a
generative-chemistry paper is expected to report. The first three are counts
and can be recomputed by anyone from the records already deposited. Internal
diversity cannot: it is a mean over pairs, the pair count grows with the square
of the set, and the largest family has tens of thousands of products, so it is
measured over a sample. A sample that lives only in the figure code is a number
a reader has to take on trust.

So it is measured here and deposited: the value, the settings that produced it,
and the identifier of every product the sample drew. Anyone can recompute the
figure from the sample, or redraw the sample from the seed and get the same one.

Each quantity is defined narrowly, because all four names are used loosely:

Validity is stronger here than the usual sense. Elsewhere it asks whether a
proposal is a molecule. Here it asks whether a proposal is a molecule *and* a
member of the family it was edited from, which is what certification tests.

Uniqueness is the share of certified products that are distinct from one
another. A product two operators both reached counts once.

Novelty is the share of distinct products that are not among the curated
ions the engine was given. It is measured against that inventory and nothing
else: it says nothing about the literature or any compound collection.

Internal diversity is one minus the mean pairwise Tanimoto similarity among
distinct products, over the fingerprint the run itself recorded.
"""
from __future__ import annotations

import csv
import json
import os
import random
import statistics
from typing import Any, Dict, List, Sequence, Tuple

from rdkit import DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from .common import RunError, parse, write_csv, write_json

RDLogger.DisableLog("rdApp.*")

#: how many distinct products of a family enter the pairwise mean
SAMPLE = 900

#: so the sample is the same every run
SEED = 20261005

COLUMNS: Sequence[str] = ("family_id", "examined", "certified", "distinct",
                          "novel", "validity", "uniqueness", "novelty",
                          "diversity", "diversity_sampled",
                          "median_similarity_to_parent")

SAMPLE_COLUMNS: Sequence[str] = ("family_id", "product_smiles")


def _settings(root: str) -> Tuple[int, int, str, int]:
    """The similarity the run measured with, and the budget it reports on."""
    path = os.path.join(root, "05_single_edits", "single_edit_settings.json")
    with open(path, encoding="utf-8") as handle:
        spec = json.load(handle)
    similarity = spec["similarity"]
    return (int(similarity["radius"]), int(similarity["bits"]),
            similarity["kind"], int(spec["allowance_the_study_reports"]))


def _given(root: str) -> set:
    """The curated ions, which is what novelty is measured against."""
    path = os.path.join(root, "02_curation", "curation_audit.csv")
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["canonical_smiles"] for row in csv.DictReader(handle)
                if row["status"] == "retained" and row["canonical_smiles"]}


def _mean_pairwise(products: Sequence[str], maker) -> float:
    """One minus the mean pairwise similarity. None when there is no pair."""
    if len(products) < 2:
        return None
    prints = []
    for smiles in products:
        mol = parse(smiles)
        if mol is None:
            raise RunError("a deposited product does not parse: %s" % smiles)
        prints.append(maker.GetFingerprint(mol))
    total, pairs = 0.0, 0
    for index in range(1, len(prints)):
        values = DataStructs.BulkTanimotoSimilarity(prints[index],
                                                    prints[:index])
        total += sum(values)
        pairs += len(values)
    return 1.0 - total / pairs if pairs else None


def _drawn(products: Sequence[str]) -> List[str]:
    """The sample the mean is taken over, drawn the same way every run."""
    ordered = sorted(products)
    if len(ordered) <= SAMPLE:
        return ordered
    return sorted(random.Random(SEED).sample(ordered, SAMPLE))


def build(root: str) -> Dict[str, Any]:
    """Measure the four quantities, deposit them and the sample they used."""
    radius, bits, kind, allowance = _settings(root)
    given = _given(root)

    path = os.path.join(root, "05_single_edits", "examination_records.csv")
    if not os.path.exists(path):
        raise RunError("the single-edit records are not in the archive")

    seen: Dict[str, set] = {}
    examined: Dict[str, int] = {}
    certified: Dict[str, int] = {}
    closeness: Dict[str, List[float]] = {}
    csv.field_size_limit(1 << 30)
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["order_index"] or 0) > allowance:
                continue
            family = row["family_id"]
            examined[family] = examined.get(family, 0) + 1
            if row["certified"] != "1" or not row["product_smiles"]:
                continue
            certified[family] = certified.get(family, 0) + 1
            seen.setdefault(family, set()).add(row["product_smiles"])
            if row.get("similarity_to_parent"):
                closeness.setdefault(family, []).append(
                    float(row["similarity_to_parent"]))

    maker = rdFingerprintGenerator.GetMorganGenerator(radius=radius,
                                                      fpSize=bits)
    rows: List[Dict[str, Any]] = []
    sampled: List[Dict[str, str]] = []

    def measure(family: str, products: Sequence[str]) -> Dict[str, Any]:
        drawn = _drawn(products)
        for smiles in drawn:
            sampled.append({"family_id": family, "product_smiles": smiles})
        value = _mean_pairwise(drawn, maker)
        fresh = sum(1 for s in products if s not in given)
        look = examined.get(family, 0) if family != "all families" else \
            sum(examined.values())
        good = certified.get(family, 0) if family != "all families" else \
            sum(certified.values())
        near = (closeness.get(family, []) if family != "all families"
                else [v for group in closeness.values() for v in group])
        middle = statistics.median(near) if near else None
        return {
            "family_id": family,
            "examined": look,
            "certified": good,
            "distinct": len(products),
            "novel": fresh,
            "validity": round(100.0 * good / look, 4) if look else None,
            "uniqueness": (round(100.0 * len(products) / good, 4)
                           if good else None),
            "novelty": (round(100.0 * fresh / len(products), 4)
                        if products else None),
            "diversity": None if value is None else round(value, 6),
            "diversity_sampled": len(drawn),
            "median_similarity_to_parent": (round(middle, 6)
                                            if middle is not None
                                            else None),
        }

    for family in sorted(seen):
        rows.append(measure(family, sorted(seen[family])))
    everything = sorted({s for group in seen.values() for s in group})
    whole = measure("all families", everything)
    rows.append(whole)

    write_csv(os.path.join(root, "10_aggregates", "product_quality.csv"),
              rows, COLUMNS)
    write_csv(os.path.join(root, "10_aggregates", "diversity_sample.csv"),
              sampled, SAMPLE_COLUMNS)
    write_json(os.path.join(root, "10_aggregates", "quality_settings.json"), {
        "allowance_measured_at": allowance,
        "similarity": {"kind": kind, "radius": radius, "bits": bits},
        "diversity_sample_ceiling_per_family": SAMPLE,
        "seed": SEED,
        "sample_rule": "the distinct products of a family, sorted, and if more "
                       "than the ceiling a seeded sample of that many; the "
                       "drawn set is deposited beside this file",
        "definitions": {
            "validity": "certified products as a share of products examined. "
                        "Certified means a molecule and a member of the family "
                        "it was edited from, which is stronger than the usual "
                        "sense of a molecule alone",
            "uniqueness": "distinct canonical products as a share of certified "
                          "products; a product two operators both reached "
                          "counts once",
            "novelty": "distinct products that are not among the curated ions "
                       "the engine was given, as a share of distinct products. "
                       "Measured against that inventory and nothing else",
            "diversity": "one minus the mean pairwise Tanimoto similarity "
                         "among the sampled distinct products",
        }})
    return {"rows": rows, "whole": whole, "sampled": len(sampled)}
