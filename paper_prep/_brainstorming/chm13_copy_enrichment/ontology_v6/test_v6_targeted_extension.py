#!/usr/bin/env python3
"""Small exact tests for selective physical-CN extension helpers."""

from __future__ import annotations

import unittest

import numpy as np

import v6_engine as V
import v6_targeted_extension as X


class TargetedExtensionTests(unittest.TestCase):
    def test_reduced_arrays_retain_cn_and_only_selected_hypotheses(self) -> None:
        full = V.WeightedTermArrays([
            np.asarray([0, 1, 4], dtype=np.int32),
            np.asarray([1, 3], dtype=np.int32),
        ], np.asarray([65, 2], dtype=np.uint32))
        reduced = X.reduced_arrays(full, [1, 4])
        np.testing.assert_array_equal(reduced[0], [0, 1])
        np.testing.assert_array_equal(reduced[1], [0])
        np.testing.assert_array_equal(reduced.weights, [65, 2])
        np.testing.assert_array_equal(V.recount(np.asarray([0, 1]), reduced, 2), [67, 65])

    def test_update_refines_raw_but_never_selected_subset_maxT(self) -> None:
        row = {
            "assignment": "midpoint", "collection": "GO_BP", "relation": "direct",
            "observed_physical_copy_burden": 2, "null_mean": 0, "null_median": 0,
            "null_q025": 0, "null_q975": 0, "null_max": 0, "count_difference": 2,
            "enrichment_ratio": 5, "raw_exceedances": 0, "raw_permutations": 3,
            "p_empirical": 0.25, "p_mc95_lower": 0, "p_mc95_upper": 1,
            "p_sequential_lower": 0, "p_sequential_upper": 1, "bh_q": 0.25,
            "by_q": 0.25, "bh_sequential_lower": 0, "bh_sequential_upper": 1,
            "collection_maxT_sequential_lower": 0.1,
            "collection_maxT_sequential_upper": 0.2,
            "global_maxT_sequential_lower": 0.1, "global_maxT_sequential_upper": 0.2,
            "collection_maxT_p": 0.15, "global_maxT_p": 0.15,
            "non_informative_constant": 0, "mc_status": "MC_UNRESOLVED",
        }
        before = (row["collection_maxT_p"], row["global_maxT_p"])
        V.COLLECTIONS = ("GO_BP",)
        try:
            X.update_candidate_statistics([row], {0: np.asarray([0, 1, 2, 0])}, 249_999)
        finally:
            V.COLLECTIONS = ("GO_BP", "GO_MF", "GO_CC", "Reactome")
        self.assertEqual((row["collection_maxT_p"], row["global_maxT_p"]), before)
        self.assertEqual(row["raw_permutations"], 4)


if __name__ == "__main__":
    unittest.main()
