#!/usr/bin/env python3
"""Build the immutable CHM13 PHR physical-gene-copy input tables.

The annotation universe is exactly the set of ``gene`` feature rows in the
specified CHM13 GFF3.  In particular, this program never collapses repeated
gene names and never creates loci on chromosomes where a GFF3 row is absent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


EXPECTED_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
EXPECTED_PHR_COUNT = 37
EXPECTED_TARGET_BP = 6_014_981
EXPECTED_GENE_COUNT = 61_312
EXPECTED_MIDPOINT_COUNT = 402
EXPECTED_OVERLAP_COUNT = 412
BUILDER_VERSION = "1.0"

GENE_COLUMNS = (
    "locus_id",
    "gene_name",
    "chromosome",
    "arm",
    "start0",
    "end0",
    "start1",
    "end1",
    "length_bp",
    "midpoint0",
    "midpoint1",
    "strand",
    "gene_biotype",
    "gff_source",
    "gff_line",
    "copy_num_id",
    "extra_copy_number",
    "in_phr_midpoint",
    "in_phr_any_overlap",
    "primary_phr_id",
    "overlap_phr_id",
    "overlap_bp",
    "record_origin",
)


@dataclass(frozen=True)
class ArmCoordinates:
    chromosome: str
    q_start0: int
    chromosome_end0: int


@dataclass(frozen=True)
class PhrInterval:
    phr_id: str
    chromosome: str
    arm: str
    start0: int
    end0: int
    sharing_chromosomes: str
    source_line: int

    @property
    def length_bp(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class GeneLocus:
    locus_id: str
    gene_name: str
    chromosome: str
    arm: str
    start0: int
    end0: int
    strand: str
    gene_biotype: str
    gff_source: str
    gff_line: int
    copy_num_id: str
    extra_copy_number: str
    midpoint_phr_id: str
    overlap_phr_id: str
    overlap_bp: int

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def end1(self) -> int:
        return self.end0

    @property
    def length_bp(self) -> int:
        return self.end0 - self.start0

    @property
    def midpoint0(self) -> int:
        # For even-length features, choose the right of the two central bases.
        return (self.start0 + self.end0) // 2

    @property
    def midpoint1(self) -> int:
        return self.midpoint0 + 1

    def as_row(self) -> tuple[object, ...]:
        return (
            self.locus_id,
            self.gene_name,
            self.chromosome,
            self.arm,
            self.start0,
            self.end0,
            self.start1,
            self.end1,
            self.length_bp,
            self.midpoint0,
            self.midpoint1,
            self.strand,
            self.gene_biotype,
            self.gff_source,
            self.gff_line,
            self.copy_num_id,
            self.extra_copy_number,
            int(bool(self.midpoint_phr_id)),
            int(bool(self.overlap_phr_id)),
            self.midpoint_phr_id or ".",
            self.overlap_phr_id or ".",
            self.overlap_bp,
            "gff3_gene_row",
        )


def chromosome_key(chromosome: str) -> int:
    try:
        return EXPECTED_CHROMS.index(chromosome)
    except ValueError as exc:
        raise ValueError(f"non-canonical CHM13 chromosome: {chromosome}") from exc


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"malformed GFF3 attribute: {item!r}")
        key, value = item.split("=", 1)
        if key in attributes:
            raise ValueError(f"duplicate GFF3 attribute key: {key!r}")
        attributes[key] = value
    return attributes


def load_arm_coordinates(
    path: Path, expected_chromosomes: Sequence[str] = EXPECTED_CHROMS
) -> dict[str, ArmCoordinates]:
    """Read p/q boundaries from a contiguous CHM13v2 cytoband BED."""
    bands: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected >=4 columns")
            chromosome, start_text, end_text, band = fields[:4]
            chromosome_key(chromosome)
            start0, end0 = int(start_text), int(end_text)
            if start0 < 0 or end0 <= start0 or band[:1] not in {"p", "q"}:
                raise ValueError(f"{path}:{line_number}: invalid cytoband row")
            bands[chromosome].append((start0, end0, band))

    if set(bands) != set(expected_chromosomes):
        raise ValueError(
            f"cytoband chromosomes differ from canonical CHM13 set: {sorted(bands)}"
        )

    result: dict[str, ArmCoordinates] = {}
    for chromosome in expected_chromosomes:
        rows = bands[chromosome]
        if rows[0][0] != 0:
            raise ValueError(f"{chromosome}: cytobands do not begin at zero")
        for previous, current in zip(rows, rows[1:]):
            if previous[1] != current[0]:
                raise ValueError(f"{chromosome}: cytobands are not contiguous")
        q_starts = [start for start, _end, band in rows if band.startswith("q")]
        if not q_starts:
            raise ValueError(f"{chromosome}: no q-arm cytoband")
        q_start0 = min(q_starts)
        if any(band.startswith("q") for _s, _e, band in rows if _e <= q_start0):
            raise ValueError(f"{chromosome}: q-arm band precedes q boundary")
        if any(band.startswith("p") for start, _e, band in rows if start >= q_start0):
            raise ValueError(f"{chromosome}: p-arm band follows q boundary")
        result[chromosome] = ArmCoordinates(chromosome, q_start0, rows[-1][1])
    return result


def load_arm_coordinates_unchecked(
    path: Path, expected_chromosomes: Sequence[str]
) -> dict[str, ArmCoordinates]:
    """Fixture-facing form with an explicitly restricted chromosome set."""
    return load_arm_coordinates(path, expected_chromosomes)


def assign_arm(chromosome: str, position0: int, arms: Mapping[str, ArmCoordinates]) -> str:
    coordinates = arms[chromosome]
    if not 0 <= position0 < coordinates.chromosome_end0:
        raise ValueError(f"{chromosome}:{position0}: position outside cytoband extent")
    suffix = "p" if position0 < coordinates.q_start0 else "q"
    return f"{chromosome}_{suffix}"


def load_phrs(
    path: Path,
    arms: Mapping[str, ArmCoordinates],
    expected_count: int = EXPECTED_PHR_COUNT,
    expected_target_bp: int | None = EXPECTED_TARGET_BP,
    expected_chromosomes: Sequence[str] | None = EXPECTED_CHROMS[:-1],
) -> list[PhrInterval]:
    intervals: list[PhrInterval] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_number}: expected exactly four columns")
            chromosome, start_text, end_text, sharing = fields
            chromosome_key(chromosome)
            sharing_names = sharing.split(",")
            if not sharing_names or any(not name for name in sharing_names):
                raise ValueError(f"{path}:{line_number}: empty sharing chromosome")
            for name in sharing_names:
                chromosome_key(name)
            start0, end0 = int(start_text), int(end_text)
            if start0 < 0 or end0 <= start0:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            if end0 > arms[chromosome].chromosome_end0:
                raise ValueError(f"{path}:{line_number}: interval exceeds chromosome")
            arm = assign_arm(chromosome, (start0 + end0) // 2, arms)
            phr_id = f"CHM13_{chromosome}_{arm[-1]}_{start0}_{end0}"
            intervals.append(
                PhrInterval(phr_id, chromosome, arm, start0, end0, sharing, line_number)
            )

    if len(intervals) != expected_count:
        raise ValueError(f"expected {expected_count} PHRs, observed {len(intervals)}")
    if len({row.phr_id for row in intervals}) != len(intervals):
        raise ValueError("duplicate PHR identifiers")
    if len({row.arm for row in intervals}) != len(intervals):
        raise ValueError("the target contains more than one PHR on an arm")
    target_bp = sum(row.length_bp for row in intervals)
    if expected_target_bp is not None and target_bp != expected_target_bp:
        raise ValueError(f"expected {expected_target_bp} target bp, observed {target_bp}")
    if expected_chromosomes is not None:
        observed_chromosomes = {row.chromosome for row in intervals}
        if observed_chromosomes != set(expected_chromosomes):
            raise ValueError(
                "target chromosomes differ from expected CHM13 subset: "
                f"{sorted(observed_chromosomes, key=chromosome_key)}"
            )
    by_chromosome: dict[str, list[PhrInterval]] = defaultdict(list)
    for row in intervals:
        by_chromosome[row.chromosome].append(row)
    for chromosome, rows in by_chromosome.items():
        rows.sort(key=lambda row: row.start0)
        for left, right in zip(rows, rows[1:]):
            if left.end0 > right.start0:
                raise ValueError(f"overlapping target PHRs on {chromosome}")
    return intervals


def load_phrs_unchecked(
    path: Path, arms: Mapping[str, ArmCoordinates], expected_count: int
) -> list[PhrInterval]:
    """Fixture-facing form with an explicitly supplied interval count."""
    return load_phrs(path, arms, expected_count, None, None)


def interval_assignments(
    chromosome: str,
    start0: int,
    end0: int,
    midpoint0: int,
    phrs_by_chromosome: Mapping[str, Sequence[PhrInterval]],
) -> tuple[str, str, int]:
    midpoint_hits = [
        phr
        for phr in phrs_by_chromosome.get(chromosome, ())
        if phr.start0 <= midpoint0 < phr.end0
    ]
    overlap_hits = [
        phr
        for phr in phrs_by_chromosome.get(chromosome, ())
        if start0 < phr.end0 and end0 > phr.start0
    ]
    if len(midpoint_hits) > 1 or len(overlap_hits) > 1:
        raise ValueError(f"gene at {chromosome}:{start0}-{end0} hits multiple PHRs")
    midpoint_id = midpoint_hits[0].phr_id if midpoint_hits else ""
    overlap_id = overlap_hits[0].phr_id if overlap_hits else ""
    overlap_bp = 0
    if overlap_hits:
        phr = overlap_hits[0]
        overlap_bp = min(end0, phr.end0) - max(start0, phr.start0)
    if midpoint_id and midpoint_id != overlap_id:
        raise AssertionError("a midpoint hit must also be an overlap hit")
    return midpoint_id, overlap_id, overlap_bp


def load_gene_loci(
    path: Path,
    arms: Mapping[str, ArmCoordinates],
    phrs: Sequence[PhrInterval],
    expected_gene_count: int | None = EXPECTED_GENE_COUNT,
    expected_midpoint_count: int | None = EXPECTED_MIDPOINT_COUNT,
    expected_overlap_count: int | None = EXPECTED_OVERLAP_COUNT,
) -> list[GeneLocus]:
    phrs_by_chromosome: dict[str, list[PhrInterval]] = defaultdict(list)
    for phr in phrs:
        phrs_by_chromosome[phr.chromosome].append(phr)

    loci: list[GeneLocus] = []
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected nine GFF3 columns")
            chromosome, source, feature_type, start_text, end_text = fields[:5]
            if feature_type != "gene":
                continue
            chromosome_key(chromosome)
            start1, end1 = int(start_text), int(end_text)
            start0, end0 = start1 - 1, end1
            if start1 < 1 or end1 < start1:
                raise ValueError(f"{path}:{line_number}: invalid GFF3 coordinates")
            strand = fields[6]
            if strand not in {"+", "-", ".", "?"}:
                raise ValueError(f"{path}:{line_number}: invalid strand {strand!r}")
            attributes = parse_attributes(fields[8])
            try:
                locus_id = attributes["ID"]
                gene_name = attributes["gene_name"]
                gene_biotype = attributes["gene_biotype"]
            except KeyError as exc:
                raise ValueError(f"{path}:{line_number}: missing required attribute {exc}") from exc
            midpoint0 = (start0 + end0) // 2
            midpoint_id, overlap_id, overlap_bp = interval_assignments(
                chromosome, start0, end0, midpoint0, phrs_by_chromosome
            )
            loci.append(
                GeneLocus(
                    locus_id=locus_id,
                    gene_name=gene_name,
                    chromosome=chromosome,
                    arm=assign_arm(chromosome, midpoint0, arms),
                    start0=start0,
                    end0=end0,
                    strand=strand,
                    gene_biotype=gene_biotype,
                    gff_source=source,
                    gff_line=line_number,
                    copy_num_id=attributes.get("copy_num_ID", "."),
                    extra_copy_number=attributes.get("extra_copy_number", "."),
                    midpoint_phr_id=midpoint_id,
                    overlap_phr_id=overlap_id,
                    overlap_bp=overlap_bp,
                )
            )

    if expected_gene_count is not None and len(loci) != expected_gene_count:
        raise ValueError(f"expected {expected_gene_count} gene rows, observed {len(loci)}")
    if expected_gene_count is not None:
        observed_chromosomes = {row.chromosome for row in loci}
        if observed_chromosomes != set(EXPECTED_CHROMS):
            raise ValueError(
                "GFF3 gene chromosomes differ from canonical CHM13 set: "
                f"{sorted(observed_chromosomes, key=chromosome_key)}"
            )
    identifiers = [row.locus_id for row in loci]
    duplicates = [name for name, count in Counter(identifiers).items() if count > 1]
    if duplicates:
        raise ValueError(f"GFF3 gene IDs are not unique: {duplicates[:5]}")
    midpoint_count = sum(bool(row.midpoint_phr_id) for row in loci)
    overlap_count = sum(bool(row.overlap_phr_id) for row in loci)
    if expected_midpoint_count is not None and midpoint_count != expected_midpoint_count:
        raise ValueError(
            f"expected {expected_midpoint_count} midpoint loci, observed {midpoint_count}"
        )
    if expected_overlap_count is not None and overlap_count != expected_overlap_count:
        raise ValueError(
            f"expected {expected_overlap_count} overlapping loci, observed {overlap_count}"
        )
    return loci


def load_gene_loci_unchecked(
    path: Path,
    arms: Mapping[str, ArmCoordinates],
    phrs: Sequence[PhrInterval],
) -> list[GeneLocus]:
    """Fixture-facing form that retains structural checks but not real-data totals."""
    return load_gene_loci(path, arms, phrs, None, None, None)


def open_deterministic_gzip_text(path: Path) -> io.TextIOWrapper:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    if path.suffix == ".gz":
        handle = open_deterministic_gzip_text(path)
    else:
        handle = path.open("w", newline="")
    with handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def write_analysis_tables(
    output_dir: Path,
    phrs: Sequence[PhrInterval],
    loci: Sequence[GeneLocus],
    arms: Mapping[str, ArmCoordinates],
    inputs: Mapping[str, Path],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    expected_files = {
        "chm13_phr_intervals.tsv",
        "chm13_gene_loci.tsv.gz",
        "chm13_phr_gene_midpoint.tsv",
        "chm13_phr_gene_any_overlap.tsv",
        "chm13_arm_summary.tsv",
        "PROVENANCE.tsv",
        "MANIFEST.sha256",
    }
    unexpected = {path.name for path in output_dir.iterdir()} - expected_files
    if unexpected:
        raise ValueError(f"refusing to mix outputs with unexpected files: {sorted(unexpected)}")

    row_counts: dict[str, int] = {}
    row_counts["chm13_phr_intervals.tsv"] = write_tsv(
        output_dir / "chm13_phr_intervals.tsv",
        (
            "phr_id",
            "chromosome",
            "arm",
            "start0",
            "end0",
            "length_bp",
            "sharing_chromosomes",
            "source_bed_line",
        ),
        (
            (
                phr.phr_id,
                phr.chromosome,
                phr.arm,
                phr.start0,
                phr.end0,
                phr.length_bp,
                phr.sharing_chromosomes,
                phr.source_line,
            )
            for phr in phrs
        ),
    )
    row_counts["chm13_gene_loci.tsv.gz"] = write_tsv(
        output_dir / "chm13_gene_loci.tsv.gz", GENE_COLUMNS, (row.as_row() for row in loci)
    )
    midpoint_loci = [row for row in loci if row.midpoint_phr_id]
    overlap_loci = [row for row in loci if row.overlap_phr_id]
    row_counts["chm13_phr_gene_midpoint.tsv"] = write_tsv(
        output_dir / "chm13_phr_gene_midpoint.tsv",
        GENE_COLUMNS,
        (row.as_row() for row in midpoint_loci),
    )
    row_counts["chm13_phr_gene_any_overlap.tsv"] = write_tsv(
        output_dir / "chm13_phr_gene_any_overlap.tsv",
        GENE_COLUMNS,
        (row.as_row() for row in overlap_loci),
    )

    locus_by_arm = Counter(row.arm for row in loci)
    midpoint_by_arm = Counter(row.arm for row in midpoint_loci)
    overlap_by_arm = Counter(row.arm for row in overlap_loci)
    phrs_by_arm = Counter(row.arm for row in phrs)
    target_bp_by_arm = Counter()
    for phr in phrs:
        target_bp_by_arm[phr.arm] += phr.length_bp
    arm_rows = []
    for chromosome in EXPECTED_CHROMS:
        coordinates = arms[chromosome]
        for suffix, start0, end0 in (
            ("p", 0, coordinates.q_start0),
            ("q", coordinates.q_start0, coordinates.chromosome_end0),
        ):
            arm = f"{chromosome}_{suffix}"
            arm_rows.append(
                (
                    arm,
                    chromosome,
                    start0,
                    end0,
                    end0 - start0,
                    locus_by_arm[arm],
                    phrs_by_arm[arm],
                    target_bp_by_arm[arm],
                    midpoint_by_arm[arm],
                    overlap_by_arm[arm],
                )
            )
    row_counts["chm13_arm_summary.tsv"] = write_tsv(
        output_dir / "chm13_arm_summary.tsv",
        (
            "arm",
            "chromosome",
            "start0",
            "end0",
            "length_bp",
            "gene_locus_count",
            "target_phr_count",
            "target_bp",
            "midpoint_gene_count",
            "any_overlap_gene_count",
        ),
        arm_rows,
    )

    provenance_rows = []
    input_roles = {
        "phr_bed": "target intervals; BED 0-based half-open",
        "gene_gff3": "physical gene loci; GFF3 1-based closed",
        "cytobands": "CHM13v2 chromosome extents and p/q boundary",
    }
    repo_root = Path(__file__).resolve().parents[3]
    for key in ("phr_bed", "gene_gff3", "cytobands"):
        path = inputs[key]
        try:
            display_path = str(path.relative_to(repo_root))
        except ValueError:
            display_path = str(path)
        provenance_rows.append(
            (key, display_path, path.stat().st_size, sha256_file(path), input_roles[key])
        )
    provenance_rows.append(
        (
            "builder",
            "build_inputs.py",
            ".",
            BUILDER_VERSION,
            "gene rows only; no symbol deduplication or locus propagation",
        )
    )
    row_counts["PROVENANCE.tsv"] = write_tsv(
        output_dir / "PROVENANCE.tsv",
        ("input", "path", "bytes", "sha256_or_version", "role"),
        provenance_rows,
    )

    manifest_rows = []
    for name in sorted(row_counts):
        path = output_dir / name
        manifest_rows.append((sha256_file(path), path.stat().st_size, row_counts[name], name))
    write_tsv(
        output_dir / "MANIFEST.sha256",
        ("sha256", "bytes", "data_rows", "file"),
        manifest_rows,
    )


def build(inputs: Mapping[str, Path], output_dir: Path) -> None:
    arms = load_arm_coordinates(inputs["cytobands"])
    phrs = load_phrs(inputs["phr_bed"], arms)
    loci = load_gene_loci(inputs["gene_gff3"], arms, phrs)
    write_analysis_tables(output_dir, phrs, loci, arms, inputs)


def compare_directories(observed: Path, rebuilt: Path) -> None:
    observed_files = sorted(path.name for path in observed.iterdir() if path.is_file())
    rebuilt_files = sorted(path.name for path in rebuilt.iterdir() if path.is_file())
    if observed_files != rebuilt_files:
        raise RuntimeError(
            f"committed/rebuilt file lists differ: {observed_files} != {rebuilt_files}"
        )
    differences = [
        name
        for name in observed_files
        if sha256_file(observed / name) != sha256_file(rebuilt / name)
    ]
    if differences:
        raise RuntimeError(f"committed outputs are not exactly reproducible: {differences}")


def default_paths() -> tuple[Path, dict[str, Path], Path]:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[2]
    inputs = {
        "phr_bed": repo_root / "data/chm13.phrs.bed",
        "gene_gff3": repo_root / "data/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz",
        "cytobands": repo_root / "data/chm13v2.0_cytobands_allchrs.bed",
    }
    return repo_root, inputs, script_dir / "analysis_ready"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    _repo_root, defaults, default_output = default_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phr-bed", type=Path, default=defaults["phr_bed"])
    parser.add_argument("--gene-gff3", type=Path, default=defaults["gene_gff3"])
    parser.add_argument("--cytobands", type=Path, default=defaults["cytobands"])
    parser.add_argument("--output-dir", type=Path, default=default_output)
    parser.add_argument(
        "--check",
        action="store_true",
        help="rebuild in a temporary directory and byte-compare with --output-dir",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = {
        "phr_bed": args.phr_bed.resolve(),
        "gene_gff3": args.gene_gff3.resolve(),
        "cytobands": args.cytobands.resolve(),
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.check:
        if not args.output_dir.is_dir():
            raise FileNotFoundError(args.output_dir)
        with tempfile.TemporaryDirectory(prefix="chm13-copy-audit-") as temporary:
            rebuilt = Path(temporary) / "analysis_ready"
            build(inputs, rebuilt)
            compare_directories(args.output_dir, rebuilt)
        print(f"PASS: {args.output_dir} is byte-for-byte reproducible")
    else:
        temporary_parent = Path(
            tempfile.mkdtemp(prefix=".analysis_ready-", dir=args.output_dir.parent)
        )
        rebuilt = temporary_parent / "analysis_ready"
        try:
            build(inputs, rebuilt)
            if args.output_dir.exists():
                shutil.rmtree(args.output_dir)
            rebuilt.replace(args.output_dir)
        finally:
            shutil.rmtree(temporary_parent, ignore_errors=True)
        print(f"Wrote audited tables to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
