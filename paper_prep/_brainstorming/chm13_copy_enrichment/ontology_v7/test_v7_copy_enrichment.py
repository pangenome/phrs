#!/usr/bin/env python3
"""Unit and release-level tests for the V7 finite-population analysis."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import math
from pathlib import Path
import unittest


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("v7_copy_enrichment", HERE / "v7_copy_enrichment.py")
assert SPEC and SPEC.loader
v7 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(v7)


class ExactInferenceTests(unittest.TestCase):
    def test_exact_upper_tail_matches_enumerated_probability(self) -> None:
        # C(3,2) C(1,0) / C(4,2) = 1/2.
        self.assertAlmostEqual(v7.exact_upper_tail(4, 3, 2, 2), 0.5, places=15)
        self.assertEqual(v7.exact_upper_tail(4, 3, 2, 0), 1.0)
        self.assertEqual(v7.exact_upper_tail(10, 0, 3, 1), 0.0)

    def test_multiplicity_adjustments_cover_complete_family(self) -> None:
        p = [0.01, 0.04, 0.2, 1.0]
        bh = v7.adjust_bh(p)
        by = v7.adjust_by(p)
        holm = v7.adjust_holm(p)
        bonf = v7.adjust_bonferroni(p)
        self.assertEqual(len(bh), len(p))
        self.assertEqual(len(by), len(p))
        self.assertEqual(len(holm), len(p))
        self.assertEqual(len(bonf), len(p))
        self.assertAlmostEqual(bh[0], 0.04)
        self.assertTrue(all(x >= y for x, y in zip(by, bh)))
        self.assertEqual(bonf, [0.04, 0.16, 0.8, 1.0])
        self.assertEqual(holm, [0.04, 0.12, 0.4, 1.0])

    def test_coordinate_copies_partition_and_are_not_source_collapsed(self) -> None:
        copies = [
            self.copy("A", "DUX4L1", "S", 1, 1),
            self.copy("B", "DUX4L2", "S", 1, 1),
            self.copy("C", "GENE3", "X", 0, 1),
            self.copy("D", "GENE4", "Y", 0, 0),
        ]
        hypotheses = [
            self.hypothesis(0, "GO_BP", "direct", "T1", 3),
            self.hypothesis(1, "GO_BP", "ancestor", "T2", 1),
        ]
        edges = [
            self.edge("A", "GO_BP", "direct", "T1"),
            self.edge("B", "GO_BP", "direct", "T1"),
            self.edge("C", "GO_BP", "direct", "T1"),
            self.edge("D", "GO_BP", "ancestor", "T2"),
        ]
        rows, audit = v7.infer_term_results(copies, edges, hypotheses)
        by_key = {(r["assignment"], r["relation"], r["term_id"]): r for r in rows}
        primary = by_key[("midpoint", "direct", "T1")]
        zero_phr = by_key[("midpoint", "ancestor", "T2")]

        self.assertEqual((primary["a_T"], primary["b_T"], primary["c_T"], primary["d_T"]),
                         (2, 0, 1, 1))
        self.assertEqual(primary["K_T_phr_burden"], 2)
        self.assertEqual(primary["M_T_genome_burden"], 3)
        self.assertAlmostEqual(primary["p_exact_upper"], 0.5)
        self.assertEqual(zero_phr["a_T"], 0)
        self.assertEqual(zero_phr["p_exact_upper"], 1.0)
        self.assertEqual(len(rows), 2 * len(hypotheses))
        self.assertEqual(audit["midpoint"]["N"], 4)
        self.assertEqual(audit["midpoint"]["K"], 2)
        self.assertEqual(audit["midpoint"]["intersection_cn"], 0)
        self.assertEqual(audit["midpoint"]["partition_sum_cn"], 4)
        # A and B share one source but contribute two coordinate-copy units.
        self.assertEqual(primary["a_T"], 2)

    def test_effect_sizes_have_explicit_zero_states(self) -> None:
        effects = v7.effect_sizes(2, 0, 0, 8)
        self.assertEqual(effects["fold_enrichment"], "inf")
        self.assertEqual(effects["odds_ratio"], "inf")
        self.assertTrue(math.isfinite(float(effects["fold_enrichment_ha"])))
        self.assertTrue(math.isfinite(float(effects["odds_ratio_ha"])))

    @staticmethod
    def copy(copy_id: str, name: str, source: str, midpoint: int, overlap: int) -> dict[str, str]:
        return {
            "copy_id": copy_id,
            "seqid": "chr1",
            "start0": str(ord(copy_id) * 10),
            "end0": str(ord(copy_id) * 10 + 5),
            "gene_name": name,
            "gene_synonyms": "",
            "functional_source_id": source,
            "functional_source_symbol": source,
            "physical_copy_cn": "1",
            "closure_term_count": "1",
            "direct_term_count": "1",
            "phr_midpoint_cn": str(midpoint),
            "phr_any_overlap_cn": str(overlap),
        }

    @staticmethod
    def hypothesis(index: int, collection: str, relation: str, term: str, burden: int) -> dict[str, str]:
        return {
            "hypothesis_index": str(index),
            "collection": collection,
            "relation": relation,
            "ontology": "GO",
            "namespace": "biological_process",
            "term_id": term,
            "term_name": term,
            "genome_physical_copy_burden": str(burden),
            "genome_arm_count": "1",
            "genome_source_count": "1",
            "multiplicity_family_id": collection,
            "hypothesis_status": "FROZEN_TESTED_NO_TARGET_FILTER",
        }

    @staticmethod
    def edge(copy_id: str, collection: str, relation: str, term: str) -> dict[str, str]:
        return {
            "copy_id": copy_id,
            "collection": collection,
            "relation": relation,
            "ontology": "GO",
            "namespace": "biological_process",
            "term_id": term,
            "physical_copy_cn": "1",
            "phr_midpoint_cn": "0",
            "phr_any_overlap_cn": "0",
        }


class FrozenReleaseTests(unittest.TestCase):
    @unittest.skipUnless((HERE / "TERM_RESULTS.tsv.gz").is_file(), "V7 release not generated yet")
    def test_generated_release_passes_validation(self) -> None:
        validation = v7.validate_release(HERE)
        self.assertEqual(validation["status"], "PASS", validation.get("errors"))
        with gzip.open(HERE / "TERM_RESULTS.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2 * 31235)
        self.assertTrue(any(row["assignment"] == "midpoint" and row["a_T"] == "0" for row in rows))


if __name__ == "__main__":
    unittest.main()
