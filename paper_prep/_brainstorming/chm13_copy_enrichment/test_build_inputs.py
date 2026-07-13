#!/usr/bin/env python3
"""Unit and committed-artifact tests for build_inputs.py."""

from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import build_inputs


class FixtureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="chm13-prep-test-")
        self.root = Path(self.temporary.name)
        self.cytobands = self.root / "cytobands.bed"
        self.bed = self.root / "target.bed"
        self.gff = self.root / "genes.gff3.gz"
        self.cytobands.write_text("chr1\t0\t500\tp11\tacen\nchr1\t500\t1000\tq11\tacen\n")
        self.bed.write_text("chr1\t100\t200\tchr1,chr2\n")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write_gff(self, rows: list[tuple[object, ...]]) -> None:
        with gzip.open(self.gff, "wt") as handle:
            handle.write("##gff-version 3\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")

    @staticmethod
    def row(
        locus_id: str,
        gene_name: str,
        start1: int,
        end1: int,
        strand: str = "+",
        biotype: str = "protein_coding",
    ) -> tuple[object, ...]:
        attrs = f"ID={locus_id};gene_name={gene_name};gene_biotype={biotype}"
        return ("chr1", "fixture", "gene", start1, end1, ".", strand, ".", attrs)

    def load_fixture(self, rows: list[tuple[object, ...]]) -> list[build_inputs.GeneLocus]:
        self.write_gff(rows)
        arms = build_inputs.load_arm_coordinates_unchecked(self.cytobands, ("chr1",))
        phrs = build_inputs.load_phrs_unchecked(self.bed, arms, expected_count=1)
        return build_inputs.load_gene_loci_unchecked(self.gff, arms, phrs)

    def test_boundary_and_coordinate_conventions(self) -> None:
        loci = self.load_fixture(
            [
                self.row("touch_left", "LEFT", 100, 100),
                self.row("first_base", "FIRST", 101, 101),
                self.row("last_base", "LAST", 200, 200),
                self.row("touch_right", "RIGHT", 201, 201),
                self.row("overlap_only", "SPAN", 91, 101),
            ]
        )
        by_id = {row.locus_id: row for row in loci}
        self.assertEqual((by_id["first_base"].start0, by_id["first_base"].end0), (100, 101))
        self.assertFalse(by_id["touch_left"].overlap_phr_id)
        self.assertTrue(by_id["first_base"].midpoint_phr_id)
        self.assertTrue(by_id["last_base"].midpoint_phr_id)
        self.assertFalse(by_id["touch_right"].overlap_phr_id)
        self.assertFalse(by_id["overlap_only"].midpoint_phr_id)
        self.assertEqual(by_id["overlap_only"].overlap_bp, 1)

    def test_strand_does_not_change_coordinates_or_membership(self) -> None:
        loci = self.load_fixture(
            [self.row("plus", "STRAND", 121, 140, "+"), self.row("minus", "STRAND", 121, 140, "-")]
        )
        plus, minus = loci
        self.assertEqual((plus.start0, plus.end0, plus.midpoint0), (minus.start0, minus.end0, minus.midpoint0))
        self.assertTrue(plus.midpoint_phr_id)
        self.assertTrue(minus.midpoint_phr_id)

    def test_multi_copy_symbol_rows_remain_independent(self) -> None:
        loci = self.load_fixture(
            [self.row("COPY_1", "COPY", 121, 130), self.row("COPY_2", "COPY", 151, 160)]
        )
        self.assertEqual([row.gene_name for row in loci], ["COPY", "COPY"])
        self.assertEqual({row.locus_id for row in loci}, {"COPY_1", "COPY_2"})

    def test_pseudogene_and_ncrna_are_retained(self) -> None:
        loci = self.load_fixture(
            [
                self.row("PSEUDO", "PSEUDO", 121, 130, biotype="pseudogene"),
                self.row("NCRNA", "NCRNA", 141, 150, biotype="lncRNA"),
            ]
        )
        self.assertEqual({row.gene_biotype for row in loci}, {"pseudogene", "lncRNA"})
        self.assertTrue(all(row.midpoint_phr_id for row in loci))

    def test_arm_assignment_uses_q_boundary(self) -> None:
        loci = self.load_fixture(
            [self.row("P_ARM", "ARM", 490, 490), self.row("Q_ARM", "ARM", 501, 501)]
        )
        self.assertEqual([row.arm for row in loci], ["chr1_p", "chr1_q"])

    def test_no_synthetic_or_propagated_loci(self) -> None:
        rows = [self.row("ONLY_A", "DUP", 121, 130), self.row("ONLY_B", "DUP", 501, 510)]
        loci = self.load_fixture(rows)
        self.assertEqual([row.locus_id for row in loci], ["ONLY_A", "ONLY_B"])
        self.assertEqual(len(loci), len(rows))
        self.assertTrue(all(row.as_row()[-1] == "gff3_gene_row" for row in loci))


class CommittedArtifactTests(unittest.TestCase):
    def test_real_inputs_reproduce_committed_tables(self) -> None:
        script = Path(__file__).resolve().parent / "build_inputs.py"
        result = subprocess.run(
            [sys.executable, str(script), "--check"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_committed_counts_and_row_origins(self) -> None:
        output = Path(__file__).resolve().parent / "analysis_ready"
        expected = {
            "chm13_phr_intervals.tsv": 37,
            "chm13_gene_loci.tsv.gz": 61_312,
            "chm13_phr_gene_midpoint.tsv": 402,
            "chm13_phr_gene_any_overlap.tsv": 412,
            "chm13_arm_summary.tsv": 48,
        }
        for name, expected_rows in expected.items():
            path = output / name
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
            self.assertEqual(len(rows), expected_rows, name)
            if "gene" in name:
                self.assertTrue(all(row["record_origin"] == "gff3_gene_row" for row in rows))

    def test_real_universe_is_one_to_one_with_gff3_gene_rows(self) -> None:
        script_dir = Path(__file__).resolve().parent
        gff = script_dir.parents[2] / "data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
        expected = []
        with gzip.open(gff, "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[2] != "gene":
                    continue
                attributes = build_inputs.parse_attributes(fields[8])
                expected.append(
                    (
                        str(line_number),
                        attributes["ID"],
                        attributes["gene_name"],
                        fields[0],
                        str(int(fields[3]) - 1),
                        fields[4],
                        fields[6],
                        attributes["gene_biotype"],
                    )
                )
        observed = []
        universe = script_dir / "analysis_ready/chm13_gene_loci.tsv.gz"
        with gzip.open(universe, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                observed.append(
                    tuple(
                        row[column]
                        for column in (
                            "gff_line",
                            "locus_id",
                            "gene_name",
                            "chromosome",
                            "start0",
                            "end0",
                            "strand",
                            "gene_biotype",
                        )
                    )
                )
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
