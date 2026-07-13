#!/usr/bin/env python3
"""Exact and integration tests for COPY_engine.py."""

from __future__ import annotations

import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import COPY_engine as engine


def write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


class EngineTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.arms_path = write(
            self.root / "arms.tsv",
            "arm\tchromosome\tstart0\tend0\n"
            "chr1_p\tchr1\t0\t10\n"
            "chr1_q\tchr1\t10\t20\n"
            "chr2_p\tchr2\t0\t10\n"
            "chr2_q\tchr2\t10\t20\n",
        )
        self.intervals_path = write(
            self.root / "intervals.tsv",
            "phr_id\tchromosome\tarm\tstart0\tend0\n"
            "target_p\tchr1\tchr1_p\t0\t2\n"
            "target_q\tchr1\tchr1_q\t18\t20\n",
        )
        # DUP occurs at two physical coordinates.  The engine must retain both.
        self.loci_path = write(
            self.root / "loci.tsv",
            "locus_id\tgene_name\tchromosome\tarm\tstart0\tend0\tmidpoint0\tgene_biotype\n"
            "L1\tDUP\tchr1\tchr1_p\t0\t1\t0\tprotein_coding\n"
            "L2\tDUP\tchr1\tchr1_p\t1\t2\t1\tprotein_coding\n"
            "L3\tOTHER\tchr1\tchr1_p\t4\t6\t5\tlncRNA\n"
            "L4\tQGENE\tchr1\tchr1_q\t18\t20\t19\tprotein_coding\n"
            "L5\tEDGE\tchr1\tchr1_q\t16\t19\t17\tprotein_coding\n",
        )
        self.terms_a = write(
            self.root / "terms_a.tsv",
            "locus_id\tterm_id\tterm_name\n"
            "L1\tT_DUP\tduplicate term\n"
            "L2\tT_DUP\tduplicate term\n"
            "L2\tT_DUP\tduplicate term\n"  # duplicate map edge counts once
            "L4\tT_Q\tq term\n",
        )
        self.terms_b = write(
            self.root / "terms_b.tsv",
            "locus_id\tterm_id\tterm_name\n"
            "L1\tT_DUP\ta distinct collection with the same term ID\n",
        )
        self.empty_terms = write(
            self.root / "empty.tsv", "locus_id\tterm_id\tterm_name\n"
        )

    def tearDown(self):
        self.temporary.cleanup()

    def objects(self):
        arms = engine.load_arms(self.arms_path)
        intervals = engine.load_intervals(self.intervals_path, arms)
        loci = engine.load_loci(self.loci_path, arms)
        genome = engine.GenomeIndex(loci, arms)
        blocks = engine.build_blocks(intervals, arms)
        return arms, intervals, loci, genome, blocks

    def test_physical_duplicate_symbols_are_not_deduplicated(self):
        arms, intervals, _loci, genome, _blocks = self.objects()
        collection = engine.load_collection("A", self.terms_a, genome, 1, 1)
        stats = engine.calculate_stats(genome, [collection], intervals, "midpoint")
        self.assertEqual(stats.burden, 3)  # L1, L2, and L4; L5 midpoint is outside
        self.assertEqual(int(stats.counts["A"][collection.term_ids.index("T_DUP")]), 2)
        self.assertEqual(len({genome.loci[i].gene_name for i in stats.selected}), 2)
        self.assertEqual(len(stats.selected), 3)

    def test_exact_candidate_boundaries_and_null_multiplicity(self):
        arms, _intervals, _loci, genome, blocks = self.objects()
        p_block = next(block for block in blocks if block.source_arm == "chr1_p")
        q_block = next(block for block in blocks if block.source_arm == "chr1_q")
        p_space = engine.make_space(p_block, arms["chr1_p"], [], "primary")
        q_space = engine.make_space(q_block, arms["chr1_q"], [], "primary")
        self.assertEqual(p_space.ranges, ((0, 8),))
        self.assertEqual(q_space.ranges, ((10, 18),))
        self.assertTrue(p_space.contains(0))
        self.assertTrue(p_space.contains(8))
        self.assertTrue(q_space.contains(10))
        self.assertTrue(q_space.contains(18))
        collection = engine.load_collection("A", self.terms_a, genome, 1, 1)
        # Enumerate the exact single-p-arm null. Starts 0 and 1 cover 2 and 1
        # physical DUP loci respectively; every other start covers zero.
        exact = []
        for start in range(9):
            placed = engine.PlacedBlock(1, p_block.block_id, "chr1_p", "chr1_p",
                                        start, p_block.components)
            intervals = engine.placed_intervals([placed], arms)
            exact.append(int(engine.calculate_stats(genome, [collection], intervals,
                                                    "midpoint").counts["A"][0]))
        self.assertEqual(exact, [2, 1, 0, 0, 0, 0, 0, 0, 0])

    def test_rigid_abutting_block_preserves_widths_and_gap(self):
        arms = engine.load_arms(self.arms_path)
        rows = [engine.Interval("a", "chr1", "chr1_p", 0, 2),
                engine.Interval("b", "chr1", "chr1_p", 2, 4)]
        blocks = engine.build_blocks(rows, arms)
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0].components, (("a", 0, 2), ("b", 2, 2)))
        placed = engine.PlacedBlock(1, blocks[0].block_id, "chr1_p", "chr1_p", 5,
                                    blocks[0].components)
        self.assertEqual([(row.start, row.end) for row in placed.intervals("chr1")],
                         [(5, 7), (7, 9)])
        first, second = self.root / "one.tsv.gz", self.root / "two.tsv.gz"
        engine.save_placement_batch(first, [[placed]], arms)
        engine.save_placement_batch(second, [[placed]], arms)
        self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_mask_matching_and_boundary_reachability(self):
        arms, _intervals, _loci, _genome, blocks = self.objects()
        block = next(block for block in blocks if block.source_arm == "chr1_p")
        # Observed [0,2) has 50% mask overlap. Only starts with 40--60% overlap
        # survive the exact fraction match; for integer bases this is only 0.
        space = engine.make_space(block, arms["chr1_p"], [(0, 1)], "primary", 100)
        self.assertEqual(space.explicit_starts.tolist(), [0])

    def test_overlap_boundary_and_empty_collection(self):
        arms, intervals, _loci, genome, _blocks = self.objects()
        empty = engine.load_collection("EMPTY", self.empty_terms, genome, 1, 1)
        midpoint = engine.calculate_stats(genome, [empty], intervals, "midpoint")
        overlap = engine.calculate_stats(genome, [empty], intervals, "overlap")
        self.assertEqual(empty.n_terms, 0)
        self.assertEqual(midpoint.counts["EMPTY"].shape, (0,))
        self.assertEqual(midpoint.burden, 3)
        self.assertEqual(overlap.burden, 4)
        # L5 overlaps [18,20) but its midpoint is outside; half-open edges do
        # not introduce L3 and no physical locus is counted twice.

    def test_zero_ties_and_clopper_pearson_boundaries(self):
        exceed, p_value, degenerate = engine.empirical_p(0, np.array([0, 1, 2]), True)
        self.assertEqual((exceed, p_value), (3, 1.0))
        self.assertFalse(degenerate)
        exceed, p_value, degenerate = engine.empirical_p(2, np.array([2, 2, 2]))
        self.assertEqual((exceed, p_value, degenerate), (3, 1.0, True))
        lower, upper = engine.clopper_pearson(0, 10)
        self.assertEqual(lower, 0.0)
        self.assertAlmostEqual(upper, 1 - 0.025 ** 0.1, places=10)
        lower, upper = engine.clopper_pearson(10, 10)
        self.assertAlmostEqual(lower, 0.025 ** 0.1, places=10)
        self.assertEqual(upper, 1.0)

    def test_bh_is_separate_and_max_t_matches_exact_hand_example(self):
        np.testing.assert_allclose(engine.bh_adjust([0.01, 0.04]), [0.02, 0.04])
        # A distinct one-hypothesis collection must keep q=p, rather than being
        # adjusted together with the preceding collection.
        np.testing.assert_allclose(engine.bh_adjust([0.04]), [0.04])
        observed = np.array([2.0, 1.0])
        null = np.array([[0.0, 1.0], [1.0, 0.0], [2.0, 0.0]])
        adjusted, exceed, z_observed, maxima = engine.max_t(observed, null)
        np.testing.assert_array_equal(exceed, [2, 1])
        np.testing.assert_allclose(adjusted, [0.75, 0.5])
        self.assertAlmostEqual(z_observed[0], math.sqrt(27 / 44), places=10)
        self.assertEqual(maxima.shape, (3,))

    def test_sampler_width_arm_and_terminal_stratum_invariants(self):
        arms, _intervals, _loci, _genome, blocks = self.objects()
        sampler = engine.RegionSampler(blocks, arms, {name: [] for name in arms},
                                       "primary", min_candidates=1)
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(12)))
        for replicate in range(100):
            for placed in sampler.sample(rng, replicate):
                source = next(block for block in blocks if block.block_id == placed.block_id)
                self.assertEqual(placed.arm, source.source_arm)
                self.assertEqual(placed.end - placed.start, source.span)
                midpoint = placed.start + source.midpoint_offset
                self.assertEqual(engine.stratum_index(engine.terminal_distance(arms[placed.arm], midpoint)),
                                 source.stratum)

    def test_terminal_and_adjacent_sensitivities_keep_geometry_separate(self):
        arms, _intervals, _loci, _genome, blocks = self.objects()
        masks = {name: [] for name in arms}
        adjacent = engine.RegionSampler(blocks, arms, masks, "adjacent", min_candidates=1)
        p_block = next(block for block in blocks if block.source_arm == "chr1_p")
        q_block = next(block for block in blocks if block.source_arm == "chr1_q")
        self.assertEqual(adjacent.spaces[(p_block.block_id, "chr1_p")].ranges, ((2, 8),))
        self.assertEqual(adjacent.spaces[(q_block.block_id, "chr1_q")].ranges, ((10, 16),))
        immediate = {row.source_arm: row.start for row in adjacent.immediate_adjacent()}
        self.assertEqual(immediate, {"chr1_p": 2, "chr1_q": 16})

        terminal = engine.RegionSampler(blocks, arms, masks, "terminal", min_candidates=1)
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence(91)))
        saw_cross_arm = False
        for replicate in range(20):
            for placed in terminal.sample(rng, replicate):
                source = next(block for block in blocks if block.block_id == placed.block_id)
                self.assertEqual(engine.partition_key(placed.arm),
                                 engine.partition_key(source.source_arm))
                self.assertEqual(placed.end - placed.start, source.span)
                saw_cross_arm |= placed.arm != source.source_arm
        self.assertTrue(saw_cross_arm)

    def cli_args(self, output: Path, permutations: int, resume: bool = False,
                 term_specs=None):
        values = ["--arms", str(self.arms_path), "--intervals", str(self.intervals_path),
                  "--loci", str(self.loci_path), "--output", str(output),
                  "--permutations", str(permutations), "--batch-size", "3",
                  "--min-candidates", "1", "--min-term-loci", "1",
                  "--min-term-arms", "1", "--seed", "12345", "--allow-pilot"]
        for specification in term_specs or ["A=%s" % self.terms_a, "B=%s" % self.terms_b]:
            values.extend(["--terms", specification])
        if resume:
            values.append("--resume")
        return values

    def test_cli_collections_are_distinct_outputs_and_resume_is_byte_exact(self):
        extended = self.root / "extended"
        direct = self.root / "direct"
        self.assertEqual(engine.main(self.cli_args(extended, 4)), 0)
        first_manifest = json.loads((extended / "run_manifest.json").read_text())
        self.assertEqual(first_manifest["completed_permutations"], 4)
        self.assertEqual(engine.main(self.cli_args(extended, 8, resume=True)), 0)
        self.assertEqual(engine.main(self.cli_args(direct, 8)), 0)
        with (extended / "term_results.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        keys = {(row["collection"], row["term_id"]) for row in rows}
        self.assertIn(("A", "T_DUP"), keys)
        self.assertIn(("B", "T_DUP"), keys)
        b_rows = [row for row in rows if row["collection"] == "B"]
        for row in b_rows:
            self.assertAlmostEqual(float(row["bh_q"]), float(row["raw_p"]))
        # Resume appends new batch files and produces the same replicate stream
        # and inferential table as an uninterrupted execution.
        self.assertEqual((extended / "term_results.tsv").read_bytes(),
                         (direct / "term_results.tsv").read_bytes())
        extended_placements = b"".join(path.read_bytes() for path in sorted((extended / "batches").glob("placements*")))
        direct_placements = b"".join(path.read_bytes() for path in sorted((direct / "batches").glob("placements*")))
        # gzip headers differ by file name/batching; decoded rows must agree.
        def decoded(directory):
            rows = []
            for path in sorted((directory / "batches").glob("placements*")):
                rows.extend(list(engine.load_placement_batch(path)))
            return [[(p.replicate, p.block_id, p.arm, p.start) for p in replicate] for replicate in rows]
        self.assertEqual(decoded(extended), decoded(direct))
        self.assertNotEqual(extended_placements, b"")
        self.assertNotEqual(direct_placements, b"")

        # Force a diagnostic trigger and verify that a family removal is paired
        # to the saved placements, removes both physical DUP copies, and is
        # labeled single-driver-sensitive rather than resampling a new null.
        arms, _intervals, _loci, genome, blocks = self.objects()
        collection = engine.load_collection("A", self.terms_a, genome, 1, 1)
        diagnostic = next(dict(row) for row in rows
                          if row["collection"] == "A" and row["term_id"] == "T_DUP"
                          and row["assignment"] == "midpoint" and row["statistic"] == "composition")
        diagnostic["bh_q"] = 0.05
        diagnostic["global_maxT_p"] = ""
        families = np.asarray(["F_DUP", "F_DUP", "F_OTHER", "F_Q", "F_EDGE"], dtype=object)
        identities = np.asarray(["I_DUP", "I_DUP", "I_OTHER", "I_Q", "I_EDGE"], dtype=object)
        leave_rows, flags = engine.leave_one_diagnostics(
            extended, json.loads((extended / "run_manifest.json").read_text()),
            [diagnostic], genome, arms, blocks, {"A": collection}, families, identities)
        family_row = next(row for row in leave_rows
                          if row["grouping"] == "family" and row["group_id"] == "F_DUP")
        self.assertEqual(family_row["observed_removed_copies"], 2)
        self.assertEqual(family_row["observed_copy_count"], 0)
        self.assertEqual(family_row["permutations"], 8)
        self.assertEqual(family_row["single_driver_sensitive"], 1)
        self.assertEqual(flags[("A", "T_DUP")]["single_driver_sensitive"], 1)
        identity_row = next(row for row in leave_rows
                            if row["grouping"] == "identity" and row["group_id"] == "I_DUP")
        self.assertEqual(identity_row["observed_removed_copies"], 2)
        observed_stats = engine.calculate_stats(
            genome, [collection], engine.placed_intervals(engine.observed_placements(blocks), arms),
            "midpoint")
        summaries, driver_groups = engine.driver_rows(
            [diagnostic], observed_stats, {"A": collection}, genome, families, identities)
        self.assertEqual(summaries[0]["observed_copies"], 2)
        symbol_driver = next(row for row in driver_groups
                             if row["grouping"] == "symbol" and row["group_id"] == "DUP")
        self.assertEqual(symbol_driver["copy_count"], 2)

    def test_cli_empty_terms_and_safe_output_refusal(self):
        output = self.root / "empty_run"
        self.assertEqual(engine.main(self.cli_args(output, 3, term_specs=["EMPTY=%s" % self.empty_terms])), 0)
        header = (output / "term_results.tsv").read_text().splitlines()[0]
        self.assertIn("collection", header)
        with self.assertRaises(SystemExit) as context:
            engine.main(self.cli_args(output, 3, term_specs=["EMPTY=%s" % self.empty_terms]))
        self.assertEqual(context.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
