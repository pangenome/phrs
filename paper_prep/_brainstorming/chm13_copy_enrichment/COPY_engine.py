#!/usr/bin/env python3
"""Copy-number-aware CHM13 PHR region-set randomization.

The observation unit in this module is always a physical ``locus_id``.  Gene
symbols are carried only for diagnostics.  The null moves complete genomic
interval blocks and then recounts the fixed annotations at their coordinates.

The command line intentionally requires named term collections and an explicit
output directory, seed, and permutation count.  Existing output is accepted
only with ``--resume`` and an identical immutable run configuration.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - exercised by the CLI error path
    raise SystemExit(
        "COPY_engine.py requires NumPy >=1.19 (and PCG64DXSM support). "
        "See README.md for the reproducible guix command."
    ) from exc


SCHEMA_VERSION = "chm13-copy-enrichment-v1"
ENGINE_VERSION = "1.0.0"
DEFAULT_SEEDS = {
    "primary": 2026071301,
    "terminal": 2026071302,
    "adjacent": 2026071303,
}
STRATA = ((0, 500_000), (500_000, 1_000_000), (1_000_000, 2_000_000),
           (2_000_000, 5_000_000), (5_000_000, None))
STATISTICS = ("copy_burden", "composition", "breadth")


def canonical_arm_key(arm: str) -> Tuple[int, int]:
    match = re.fullmatch(r"chr([0-9]+|X|Y)_([pq])", arm)
    if not match:
        raise ValueError("invalid canonical arm name: %s" % arm)
    chrom, side = match.groups()
    chrom_number = 23 if chrom == "X" else 24 if chrom == "Y" else int(chrom)
    return chrom_number, 0 if side == "p" else 1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(str(path), "rt") if path.suffix == ".gz" else path.open()


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    os.replace(str(temporary), str(path))


@dataclass(frozen=True)
class Arm:
    name: str
    chromosome: str
    start: int
    end: int

    @property
    def side(self) -> str:
        return self.name[-1]

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Interval:
    interval_id: str
    chromosome: str
    arm: str
    start: int
    end: int

    @property
    def width(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class Block:
    block_id: str
    source_arm: str
    source_start: int
    source_end: int
    # (interval_id, offset from source_start, width), including overlaps.
    components: Tuple[Tuple[str, int, int], ...]
    stratum: int

    @property
    def span(self) -> int:
        return self.source_end - self.source_start

    @property
    def midpoint_offset(self) -> int:
        return self.span // 2


@dataclass(frozen=True)
class PlacedBlock:
    replicate: int
    block_id: str
    source_arm: str
    arm: str
    start: int
    components: Tuple[Tuple[str, int, int], ...]

    @property
    def end(self) -> int:
        return self.start + max(offset + width for _name, offset, width in self.components)

    def intervals(self, chromosome: str) -> Iterator[Interval]:
        for interval_id, offset, width in self.components:
            yield Interval(interval_id, chromosome, self.arm, self.start + offset,
                           self.start + offset + width)


@dataclass(frozen=True)
class Locus:
    index: int
    locus_id: str
    gene_name: str
    chromosome: str
    arm: str
    start: int
    end: int
    midpoint: int
    biotype: str


@dataclass
class Collection:
    name: str
    source: Path
    term_ids: List[str]
    term_names: List[str]
    locus_terms: List[np.ndarray]
    genome_locus_counts: np.ndarray
    genome_arm_counts: np.ndarray
    annotated: np.ndarray
    filtered_rows: List[Mapping[str, object]]

    @property
    def n_terms(self) -> int:
        return len(self.term_ids)


@dataclass
class PlacementSpace:
    block: Block
    destination_arm: Arm
    ranges: Tuple[Tuple[int, int], ...]  # inclusive integer starts
    explicit_starts: Optional[np.ndarray] = None

    @property
    def count(self) -> int:
        if self.explicit_starts is not None:
            return int(self.explicit_starts.size)
        return sum(end - start + 1 for start, end in self.ranges)

    def sample(self, rng: np.random.Generator) -> int:
        count = self.count
        if count == 0:
            raise ValueError("cannot sample an empty placement space")
        rank = int(rng.integers(0, count))
        if self.explicit_starts is not None:
            return int(self.explicit_starts[rank])
        for start, end in self.ranges:
            width = end - start + 1
            if rank < width:
                return start + rank
            rank -= width
        raise AssertionError("placement rank fell outside ranges")

    def contains(self, start: int) -> bool:
        if self.explicit_starts is not None:
            position = int(np.searchsorted(self.explicit_starts, start))
            return position < self.explicit_starts.size and int(self.explicit_starts[position]) == start
        return any(left <= start <= right for left, right in self.ranges)


@dataclass
class AssignmentStats:
    burden: int
    selected: np.ndarray
    counts: Dict[str, np.ndarray]
    denominators: Dict[str, int]
    breadth: Dict[str, np.ndarray]


def load_arms(path: Path) -> Dict[str, Arm]:
    arms: Dict[str, Arm] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"arm", "chromosome", "start0", "end0"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("arm table lacks required columns: %s" % sorted(required))
        for row in reader:
            arm = Arm(row["arm"], row["chromosome"], int(row["start0"]), int(row["end0"]))
            canonical_arm_key(arm.name)
            if arm.end <= arm.start or arm.name in arms:
                raise ValueError("invalid or duplicate arm: %s" % arm.name)
            arms[arm.name] = arm
    return dict(sorted(arms.items(), key=lambda item: canonical_arm_key(item[0])))


def load_intervals(path: Path, arms: Mapping[str, Arm]) -> List[Interval]:
    intervals: List[Interval] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"phr_id", "chromosome", "arm", "start0", "end0"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("interval table lacks required columns: %s" % sorted(required))
        for row in reader:
            interval = Interval(row["phr_id"], row["chromosome"], row["arm"],
                                int(row["start0"]), int(row["end0"]))
            arm = arms.get(interval.arm)
            if arm is None or arm.chromosome != interval.chromosome:
                raise ValueError("interval has an unknown/incompatible arm: %s" % interval.interval_id)
            if interval.start < arm.start or interval.end > arm.end or interval.end <= interval.start:
                raise ValueError("interval lies outside its arm: %s" % interval.interval_id)
            intervals.append(interval)
    if len({row.interval_id for row in intervals}) != len(intervals):
        raise ValueError("duplicate interval identifiers")
    if not intervals:
        raise ValueError("target interval table is empty")
    return intervals


def load_loci(path: Path, arms: Mapping[str, Arm]) -> List[Locus]:
    loci: List[Locus] = []
    seen: Set[str] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"locus_id", "gene_name", "chromosome", "arm", "start0", "end0", "midpoint0"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("locus table lacks required columns: %s" % sorted(required))
        for row in reader:
            locus_id = row["locus_id"]
            if locus_id in seen:
                raise ValueError("duplicate physical locus_id: %s" % locus_id)
            seen.add(locus_id)
            arm = arms.get(row["arm"])
            start, end, midpoint = int(row["start0"]), int(row["end0"]), int(row["midpoint0"])
            if arm is None or row["chromosome"] != arm.chromosome:
                raise ValueError("locus has an unknown/incompatible arm: %s" % locus_id)
            if start < arm.start or end > arm.end or not start <= midpoint < end:
                raise ValueError("locus coordinates are invalid: %s" % locus_id)
            loci.append(Locus(len(loci), locus_id, row["gene_name"], row["chromosome"],
                              row["arm"], start, end, midpoint, row.get("gene_biotype", ".")))
    if not loci:
        raise ValueError("physical locus universe is empty")
    return loci


def terminal_distance(arm: Arm, midpoint: int) -> int:
    return midpoint - arm.start if arm.side == "p" else arm.end - midpoint


def stratum_index(distance: int) -> int:
    if distance < 0:
        raise ValueError("negative terminal distance")
    for index, (lower, upper) in enumerate(STRATA):
        if distance >= lower and (upper is None or distance < upper):
            return index
    raise AssertionError("no terminal stratum")


def build_blocks(intervals: Sequence[Interval], arms: Mapping[str, Arm]) -> List[Block]:
    """Merge overlapping/abutting same-arm templates into rigid blocks."""
    by_arm: MutableMapping[str, List[Interval]] = defaultdict(list)
    for interval in intervals:
        by_arm[interval.arm].append(interval)
    blocks: List[Block] = []
    for arm_name in sorted(by_arm, key=canonical_arm_key):
        rows = sorted(by_arm[arm_name], key=lambda row: (row.start, row.end, row.interval_id))
        groups: List[List[Interval]] = []
        for row in rows:
            if not groups or row.start > max(item.end for item in groups[-1]):
                groups.append([row])
            else:
                groups[-1].append(row)
        for group_index, group in enumerate(groups, 1):
            start = min(row.start for row in group)
            end = max(row.end for row in group)
            components = tuple((row.interval_id, row.start - start, row.width) for row in group)
            arm = arms[arm_name]
            distance = terminal_distance(arm, start + (end - start) // 2)
            blocks.append(Block("%s.block%02d" % (arm_name, group_index), arm_name,
                                start, end, components, stratum_index(distance)))
    return blocks


def load_masks(path: Optional[Path], arms: Mapping[str, Arm]) -> Dict[str, List[Tuple[int, int]]]:
    masks: Dict[str, List[Tuple[int, int]]] = {name: [] for name in arms}
    if path is None:
        return masks
    chromosome_arms: Dict[str, List[Arm]] = defaultdict(list)
    for arm in arms.values():
        chromosome_arms[arm.chromosome].append(arm)
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError("%s:%d: mask BED needs three columns" % (path, line_number))
            chromosome, start_text, end_text = fields[:3]
            start, end = int(start_text), int(end_text)
            if end <= start:
                raise ValueError("%s:%d: invalid mask interval" % (path, line_number))
            if chromosome not in chromosome_arms:
                raise ValueError("%s:%d: unknown mask chromosome" % (path, line_number))
            for arm in chromosome_arms[chromosome]:
                left, right = max(start, arm.start), min(end, arm.end)
                if left < right:
                    masks[arm.name].append((left, right))
    for name, rows in masks.items():
        masks[name] = coalesce_pairs(rows)
    return masks


def coalesce_pairs(rows: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    result: List[List[int]] = []
    for start, end in sorted(rows):
        if not result or start > result[-1][1]:
            result.append([start, end])
        else:
            result[-1][1] = max(result[-1][1], end)
    return [(start, end) for start, end in result]


def block_intervals(block: Block, start: int) -> List[Tuple[int, int]]:
    return coalesce_pairs((start + offset, start + offset + width)
                          for _name, offset, width in block.components)


def overlap_length(rows: Sequence[Tuple[int, int]], masks: Sequence[Tuple[int, int]]) -> int:
    total = 0
    for start, end in rows:
        for mask_start, mask_end in masks:
            if mask_start >= end:
                break
            if mask_end > start:
                total += max(0, min(end, mask_end) - max(start, mask_start))
    return total


def base_candidate_range(block: Block, destination: Arm) -> Tuple[int, int]:
    """Inclusive start range preserving the block's fixed terminal stratum."""
    lower_distance, upper_distance = STRATA[block.stratum]
    upper_distance = destination.length if upper_distance is None else upper_distance
    offset = block.midpoint_offset
    if destination.side == "p":
        lower = destination.start + lower_distance - offset
        upper = destination.start + upper_distance - offset - 1
    else:
        lower = destination.end - offset - upper_distance + 1
        upper = destination.end - offset - lower_distance
    lower = max(lower, destination.start)
    upper = min(upper, destination.end - block.span)
    return lower, upper


def make_space(block: Block, destination: Arm, masks: Sequence[Tuple[int, int]],
               mode: str = "primary", enumerate_limit: int = 10_000_000) -> PlacementSpace:
    if mode in {"primary", "terminal"}:
        lower, upper = base_candidate_range(block, destination)
    elif mode == "adjacent":
        source_arm = destination
        if source_arm.name != block.source_arm:
            raise ValueError("adjacent placement cannot change arm")
        if source_arm.side == "p":
            lower = block.source_end
            upper = min(source_arm.end, block.source_end + 5_000_000) - block.span
        else:
            lower = max(source_arm.start, block.source_start - 5_000_000)
            upper = block.source_start - block.span
    else:
        raise ValueError("unknown placement mode: %s" % mode)
    if upper < lower:
        return PlacementSpace(block, destination, tuple())
    if not masks:
        return PlacementSpace(block, destination, ((lower, upper),))

    # The specification permits a non-empty frozen mask only when candidate and
    # template excluded fractions agree within one percentage point.  Exact
    # enumeration is deliberately used here: it is auditable and cannot bias
    # boundary starts.  Large masked spaces must be precomputed upstream.
    possible = upper - lower + 1
    if possible > enumerate_limit:
        raise ValueError(
            "masked candidate space for %s has %d starts (limit %d); precompute "
            "an exact candidate table or increase --mask-enumeration-limit"
            % (block.block_id, possible, enumerate_limit)
        )
    observed_bp = overlap_length(block_intervals(block, block.source_start), masks)
    observed_total = sum(end - start for start, end in block_intervals(block, block.source_start))
    observed_fraction = observed_bp / observed_total
    valid = []
    for start in range(lower, upper + 1):
        rows = block_intervals(block, start)
        candidate_fraction = overlap_length(rows, masks) / sum(end - left for left, end in rows)
        if abs(candidate_fraction - observed_fraction) <= 0.01 + 1e-15:
            valid.append(start)
    return PlacementSpace(block, destination, tuple(), np.asarray(valid, dtype=np.int64))


def rows_overlap(first: Sequence[Tuple[int, int]], second: Sequence[Tuple[int, int]]) -> bool:
    return any(a < d and c < b for a, b in first for c, d in second)


def partition_key(arm_name: str) -> Tuple[str, str, str]:
    chrom = arm_name[3:].split("_", 1)[0]
    acrocentric = chrom in {"13", "14", "15", "21", "22"}
    sex = chrom in {"X", "Y"}
    return arm_name[-1], "acrocentric" if acrocentric else "non_acrocentric", "sex" if sex else "autosome"


class RegionSampler:
    """Exactly uniform integer-translation sampler for joint region sets."""

    def __init__(self, blocks: Sequence[Block], arms: Mapping[str, Arm],
                 masks: Mapping[str, Sequence[Tuple[int, int]]], mode: str,
                 min_candidates: int = 100, enumeration_limit: int = 10_000_000):
        self.blocks = list(sorted(blocks, key=lambda row: (canonical_arm_key(row.source_arm), row.block_id)))
        self.arms = arms
        self.masks = masks
        self.mode = mode
        self.min_candidates = min_candidates
        self.enumeration_limit = enumeration_limit
        self.rejected_by_arm: Counter = Counter()
        self.accepted_by_arm: Counter = Counter()
        self.cross_arm_fallback: Set[Tuple[str, str, str]] = set()
        self.spaces: Dict[Tuple[str, str], PlacementSpace] = {}
        if mode in {"primary", "adjacent"}:
            for block in self.blocks:
                arm = arms[block.source_arm]
                self.spaces[(block.block_id, arm.name)] = make_space(
                    block, arm, masks.get(arm.name, ()), mode, enumeration_limit)
        elif mode == "terminal":
            for block in self.blocks:
                for arm in arms.values():
                    if partition_key(arm.name) == partition_key(block.source_arm):
                        self.spaces[(block.block_id, arm.name)] = make_space(
                            block, arm, masks.get(arm.name, ()), mode, enumeration_limit)
        else:
            raise ValueError("unknown sampler mode: %s" % mode)
        self._validate_spaces()

    def _validate_spaces(self) -> None:
        if self.mode in {"primary", "adjacent"}:
            failures = [(block.block_id, self.spaces[(block.block_id, block.source_arm)].count)
                        for block in self.blocks
                        if self.spaces[(block.block_id, block.source_arm)].count < self.min_candidates]
            if failures:
                raise ValueError("non-estimable %s candidate spaces (<%d): %s" %
                                 (self.mode, self.min_candidates, failures))
        else:
            by_partition: MutableMapping[Tuple[str, str, str], List[Block]] = defaultdict(list)
            for block in self.blocks:
                by_partition[partition_key(block.source_arm)].append(block)
            for key, blocks in by_partition.items():
                eligible = [arm for arm in self.arms.values() if partition_key(arm.name) == key]
                if len(eligible) < len(blocks):
                    self.cross_arm_fallback.add(key)
                for block in blocks:
                    if not any(self.spaces[(block.block_id, arm.name)].count >= self.min_candidates
                               for arm in eligible):
                        raise ValueError("terminal block has no estimable destination: %s" % block.block_id)

    def candidate_rows(self) -> Iterator[Mapping[str, object]]:
        for (block_id, arm_name), space in sorted(self.spaces.items()):
            yield {
                "mode": self.mode,
                "block_id": block_id,
                "source_arm": space.block.source_arm,
                "destination_arm": arm_name,
                "candidate_count": space.count,
                "range_count": len(space.ranges),
                "explicit": int(space.explicit_starts is not None),
                "observed_start_is_candidate": int(
                    arm_name == space.block.source_arm and space.contains(space.block.source_start)
                ),
            }

    def _sample_same_arm(self, rng: np.random.Generator, replicate: int) -> List[PlacedBlock]:
        by_arm: MutableMapping[str, List[Block]] = defaultdict(list)
        for block in self.blocks:
            by_arm[block.source_arm].append(block)
        result: List[PlacedBlock] = []
        for arm_name in sorted(by_arm, key=canonical_arm_key):
            blocks = by_arm[arm_name]
            while True:
                starts = [self.spaces[(block.block_id, arm_name)].sample(rng) for block in blocks]
                interval_sets = [block_intervals(block, start) for block, start in zip(blocks, starts)]
                if not any(rows_overlap(interval_sets[i], interval_sets[j])
                           for i in range(len(blocks)) for j in range(i + 1, len(blocks))):
                    break
                self.rejected_by_arm[arm_name] += 1
            self.accepted_by_arm[arm_name] += 1
            for block, start in zip(blocks, starts):
                result.append(PlacedBlock(replicate, block.block_id, block.source_arm,
                                          arm_name, start, block.components))
        return result

    def _sample_terminal(self, rng: np.random.Generator, replicate: int) -> List[PlacedBlock]:
        by_partition: MutableMapping[Tuple[str, str, str], List[Block]] = defaultdict(list)
        for block in self.blocks:
            by_partition[partition_key(block.source_arm)].append(block)
        result: List[PlacedBlock] = []
        for key in sorted(by_partition):
            blocks = sorted(by_partition[key], key=lambda row: row.block_id)
            eligible = sorted((arm for arm in self.arms.values() if partition_key(arm.name) == key),
                              key=lambda arm: canonical_arm_key(arm.name))
            fallback = key in self.cross_arm_fallback
            attempts = 0
            while True:
                attempts += 1
                if attempts > 1_000_000:
                    raise RuntimeError("terminal arm assignment acceptance fell below safe bounds: %s" % (key,))
                if fallback:
                    destinations = [self.arms[block.source_arm] for block in blocks]
                else:
                    order = rng.permutation(len(eligible))[:len(blocks)]
                    destinations = [eligible[int(index)] for index in order]
                if all(self.spaces[(block.block_id, arm.name)].count >= self.min_candidates
                       for block, arm in zip(blocks, destinations)):
                    break
            for block, arm in zip(blocks, destinations):
                space = self.spaces[(block.block_id, arm.name)]
                start = space.sample(rng)
                self.accepted_by_arm[arm.name] += 1
                result.append(PlacedBlock(replicate, block.block_id, block.source_arm,
                                          arm.name, start, block.components))
        return result

    def sample(self, rng: np.random.Generator, replicate: int) -> List[PlacedBlock]:
        if self.mode == "terminal":
            return self._sample_terminal(rng, replicate)
        return self._sample_same_arm(rng, replicate)

    def immediate_adjacent(self) -> List[PlacedBlock]:
        if self.mode != "adjacent":
            raise ValueError("immediate comparator exists only for adjacent mode")
        rows = []
        for block in self.blocks:
            space = self.spaces[(block.block_id, block.source_arm)]
            if space.count == 0:
                continue
            if space.explicit_starts is not None:
                start = int(space.explicit_starts[0] if self.arms[block.source_arm].side == "p"
                            else space.explicit_starts[-1])
            else:
                start = space.ranges[0][0] if self.arms[block.source_arm].side == "p" else space.ranges[-1][1]
            rows.append(PlacedBlock(-1, block.block_id, block.source_arm, block.source_arm,
                                    start, block.components))
        return rows


class GenomeIndex:
    def __init__(self, loci: Sequence[Locus], arms: Mapping[str, Arm]):
        self.loci = list(loci)
        self.arms = arms
        self.by_id = {locus.locus_id: locus.index for locus in loci}
        self.arm_indices: Dict[str, np.ndarray] = {}
        self.arm_midpoints: Dict[str, np.ndarray] = {}
        self.arm_starts: Dict[str, np.ndarray] = {}
        self.arm_ends: Dict[str, np.ndarray] = {}
        for arm_name in arms:
            indices = sorted((locus.index for locus in loci if locus.arm == arm_name),
                             key=lambda index: (loci[index].midpoint, loci[index].start,
                                                loci[index].end, loci[index].locus_id))
            array = np.asarray(indices, dtype=np.int64)
            self.arm_indices[arm_name] = array
            self.arm_midpoints[arm_name] = np.asarray([loci[i].midpoint for i in indices], dtype=np.int64)
            # Overlap queries need an independent start-sorted index.
            start_indices = sorted((locus.index for locus in loci if locus.arm == arm_name),
                                   key=lambda index: (loci[index].start, loci[index].end, loci[index].locus_id))
            self.arm_starts[arm_name] = np.asarray(start_indices, dtype=np.int64)
            self.arm_ends[arm_name] = np.asarray([loci[i].start for i in start_indices], dtype=np.int64)

    def select(self, intervals: Sequence[Interval], assignment: str) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        by_arm: MutableMapping[str, List[Tuple[int, int]]] = defaultdict(list)
        for interval in intervals:
            by_arm[interval.arm].append((interval.start, interval.end))
        selected_by_arm: Dict[str, np.ndarray] = {}
        all_selected: List[np.ndarray] = []
        for arm_name, raw_rows in by_arm.items():
            rows = coalesce_pairs(raw_rows)
            pieces: List[np.ndarray] = []
            if assignment == "midpoint":
                positions = self.arm_midpoints[arm_name]
                indices = self.arm_indices[arm_name]
                for start, end in rows:
                    left = int(np.searchsorted(positions, start, side="left"))
                    right = int(np.searchsorted(positions, end, side="left"))
                    pieces.append(indices[left:right])
            elif assignment == "overlap":
                indices = self.arm_starts[arm_name]
                starts = self.arm_ends[arm_name]
                for start, end in rows:
                    right = int(np.searchsorted(starts, end, side="left"))
                    candidates = indices[:right]
                    pieces.append(np.asarray([i for i in candidates if self.loci[int(i)].end > start],
                                             dtype=np.int64))
            else:
                raise ValueError("unknown locus assignment: %s" % assignment)
            combined = np.unique(np.concatenate(pieces)) if pieces else np.empty(0, dtype=np.int64)
            selected_by_arm[arm_name] = combined
            all_selected.append(combined)
        selected = np.unique(np.concatenate(all_selected)) if all_selected else np.empty(0, dtype=np.int64)
        return selected, selected_by_arm


def load_collection(name: str, path: Path, genome: GenomeIndex,
                    min_loci: int, min_arms: int) -> Collection:
    raw: MutableMapping[str, Set[int]] = defaultdict(set)
    names: Dict[str, str] = {}
    unknown_ids: Set[str] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"locus_id", "term_id"}.issubset(reader.fieldnames):
            raise ValueError("term map %s requires locus_id and term_id" % path)
        for row_number, row in enumerate(reader, 2):
            locus_id, term_id = row["locus_id"], row["term_id"]
            if not term_id:
                raise ValueError("%s:%d: empty term_id" % (path, row_number))
            if locus_id not in genome.by_id:
                unknown_ids.add(locus_id)
                continue
            raw[term_id].add(genome.by_id[locus_id])
            term_name = row.get("term_name", term_id) or term_id
            if term_id in names and names[term_id] != term_name:
                raise ValueError("term %s has conflicting names in %s" % (term_id, path))
            names[term_id] = term_name
    filtered_rows: List[Mapping[str, object]] = []
    retained: List[str] = []
    for term_id in sorted(raw):
        loci = raw[term_id]
        arms = {genome.loci[index].arm for index in loci}
        reason = "retained"
        if len(loci) < min_loci:
            reason = "genome_loci_below_%d" % min_loci
        elif len(arms) < min_arms:
            reason = "genome_arms_below_%d" % min_arms
        else:
            retained.append(term_id)
        filtered_rows.append({"collection": name, "term_id": term_id,
                              "term_name": names[term_id], "genome_loci": len(loci),
                              "genome_arms": len(arms), "status": reason})
    term_position = {term_id: index for index, term_id in enumerate(retained)}
    locus_terms: List[List[int]] = [[] for _ in genome.loci]
    genome_counts = np.zeros(len(retained), dtype=np.int64)
    genome_arm_counts = np.zeros(len(retained), dtype=np.int64)
    for term_id in retained:
        term_index = term_position[term_id]
        indices = raw[term_id]
        genome_counts[term_index] = len(indices)
        genome_arm_counts[term_index] = len({genome.loci[index].arm for index in indices})
        for locus_index in indices:
            locus_terms[locus_index].append(term_index)
    arrays = [np.asarray(sorted(values), dtype=np.int32) for values in locus_terms]
    annotated = np.asarray([values.size > 0 for values in arrays], dtype=np.bool_)
    if unknown_ids:
        examples = sorted(unknown_ids)[:5]
        raise ValueError("term map %s contains %d locus IDs outside the frozen physical universe; examples: %s"
                         % (path, len(unknown_ids), examples))
    return Collection(name, path, retained, [names[term] for term in retained], arrays,
                      genome_counts, genome_arm_counts, annotated, filtered_rows)


def load_group_map(path: Optional[Path], genome: GenomeIndex, label: str) -> np.ndarray:
    values: List[Optional[str]] = [None] * len(genome.loci)
    if path is not None:
        with open_text(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            group_column = "%s_id" % label
            if not reader.fieldnames or not {"locus_id", group_column}.issubset(reader.fieldnames):
                raise ValueError("%s map requires locus_id and %s" % (label, group_column))
            for row in reader:
                locus_id = row["locus_id"]
                if locus_id not in genome.by_id:
                    continue
                index = genome.by_id[locus_id]
                group = row[group_column]
                if not group:
                    continue
                if values[index] is not None and values[index] != group:
                    raise ValueError("locus %s has multiple %s groups" % (locus_id, label))
                values[index] = group
    # Missing identifiers are explicitly singleton groups, never one shared NA group.
    return np.asarray([value if value is not None else "__singleton__:%s" % locus.locus_id
                       for value, locus in zip(values, genome.loci)], dtype=object)


def placed_intervals(placements: Sequence[PlacedBlock], arms: Mapping[str, Arm]) -> List[Interval]:
    rows: List[Interval] = []
    for placed in placements:
        rows.extend(placed.intervals(arms[placed.arm].chromosome))
    return rows


def observed_placements(blocks: Sequence[Block]) -> List[PlacedBlock]:
    return [PlacedBlock(0, block.block_id, block.source_arm, block.source_arm,
                        block.source_start, block.components) for block in blocks]


def calculate_stats(genome: GenomeIndex, collections: Sequence[Collection],
                    intervals: Sequence[Interval], assignment: str) -> AssignmentStats:
    selected, selected_by_arm = genome.select(intervals, assignment)
    counts: Dict[str, np.ndarray] = {}
    breadth: Dict[str, np.ndarray] = {}
    denominators: Dict[str, int] = {}
    for collection in collections:
        arrays = [collection.locus_terms[int(index)] for index in selected
                  if collection.locus_terms[int(index)].size]
        annotations = np.concatenate(arrays) if arrays else np.empty(0, dtype=np.int32)
        counts[collection.name] = np.bincount(annotations, minlength=collection.n_terms).astype(np.int32)
        denominators[collection.name] = int(np.count_nonzero(collection.annotated[selected]))
        term_breadth = np.zeros(collection.n_terms, dtype=np.uint16)
        for arm_selected in selected_by_arm.values():
            arm_arrays = [collection.locus_terms[int(index)] for index in arm_selected
                          if collection.locus_terms[int(index)].size]
            if arm_arrays:
                present = np.unique(np.concatenate(arm_arrays))
                term_breadth[present] += 1
        breadth[collection.name] = term_breadth
    return AssignmentStats(int(selected.size), selected, counts, denominators, breadth)


def empirical_p(observed: float, null: np.ndarray, zero_is_one: bool = False) -> Tuple[int, float, bool]:
    null = np.asarray(null)
    degenerate = bool(null.size and np.all(null == observed))
    if zero_is_one and observed == 0:
        return int(null.size), 1.0, degenerate
    if degenerate:
        return int(null.size), 1.0, True
    exceedances = int(np.count_nonzero(null >= observed))
    return exceedances, (exceedances + 1.0) / (null.size + 1.0), False


def _betacf(a: float, b: float, x: float) -> float:
    # Modified Lentz continued fraction (Numerical Recipes), used only for
    # exact binomial interval inversion.  Double precision is ample here.
    max_iterations, epsilon, floor = 300, 3.0e-14, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = floor if abs(d) < floor else d
    d = 1.0 / d
    h = d
    for iteration in range(1, max_iterations + 1):
        m2 = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = floor if abs(d) < floor else d
        c = 1.0 + aa / c
        c = floor if abs(c) < floor else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = floor if abs(d) < floor else d
        c = 1.0 + aa / c
        c = floor if abs(c) < floor else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            return h
    raise ArithmeticError("incomplete-beta continued fraction did not converge")


def regularized_beta(x: float, a: float, b: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_factor = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    factor = math.exp(log_factor + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return factor * _betacf(a, b, x) / a
    return 1.0 - factor * _betacf(b, a, 1.0 - x) / b


def beta_quantile(probability: float, a: float, b: float) -> float:
    if probability <= 0:
        return 0.0
    if probability >= 1:
        return 1.0
    lower, upper = 0.0, 1.0
    for _ in range(100):
        middle = (lower + upper) / 2.0
        if regularized_beta(middle, a, b) < probability:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2.0


def clopper_pearson(exceedances: int, permutations: int, alpha: float = 0.05) -> Tuple[float, float]:
    if not 0 <= exceedances <= permutations:
        raise ValueError("invalid binomial count")
    lower = 0.0 if exceedances == 0 else beta_quantile(alpha / 2.0, exceedances,
                                                        permutations - exceedances + 1)
    upper = 1.0 if exceedances == permutations else beta_quantile(
        1.0 - alpha / 2.0, exceedances + 1, permutations - exceedances)
    return lower, upper


def bh_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    if values.size == 0:
        return values.copy()
    order = np.argsort(values, kind="stable")
    adjusted = np.empty(values.size, dtype=float)
    running = 1.0
    for reverse_rank in range(values.size - 1, -1, -1):
        index = int(order[reverse_rank])
        rank = reverse_rank + 1
        running = min(running, values[index] * values.size / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def by_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    harmonic = sum(1.0 / index for index in range(1, len(values) + 1)) if len(values) else 1.0
    return np.minimum(1.0, bh_adjust(values) * harmonic)


def standardized_values(observed: np.ndarray, null: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return symmetric pooled z-observed, z-null, and degeneracy flags."""
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    if null.ndim == 1:
        null = null[:, None]
    z_observed = np.zeros(observed.size, dtype=float)
    z_null = np.zeros_like(null, dtype=float)
    degenerate = np.zeros(observed.size, dtype=np.bool_)
    for index in range(observed.size):
        column = null[:, index]
        if np.all(column == observed[index]):
            degenerate[index] = True
            continue
        if np.all(column == column[0]) and observed[index] > column[0]:
            z_observed[index] = np.inf
            z_null[:, index] = 0.0
            continue
        pooled = np.concatenate(([observed[index]], column))
        standard_deviation = float(np.std(pooled, ddof=1))
        if standard_deviation == 0:
            degenerate[index] = True
            continue
        mean = float(np.mean(pooled))
        z_observed[index] = (observed[index] - mean) / standard_deviation
        z_null[:, index] = (column - mean) / standard_deviation
    return z_observed, z_null, degenerate


def max_t(observed: np.ndarray, null: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Single-step Westfall--Young maxT for one declared family."""
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    if observed.size == 0:
        return np.empty(0), np.empty(0, dtype=int), np.empty(0), np.empty(null.shape[0])
    z_observed, z_null, degenerate = standardized_values(observed, null)
    maxima = np.max(z_null, axis=1)
    exceedances = np.asarray([np.count_nonzero(maxima >= value) for value in z_observed], dtype=int)
    adjusted = (exceedances + 1.0) / (null.shape[0] + 1.0)
    adjusted[degenerate] = 1.0
    exceedances[degenerate] = null.shape[0]
    return adjusted, exceedances, z_observed, maxima


def composition_effect(observed_count: float, observed_denominator: float,
                       null_count_median: float, null_denominator_median: float) -> Tuple[float, float]:
    observed_odds = (observed_count + 0.5) / (observed_denominator - observed_count + 0.5)
    null_odds = (null_count_median + 0.5) / (null_denominator_median - null_count_median + 0.5)
    ratio = observed_odds / null_odds
    return ratio, math.log(ratio, 2)


def jsonable_rng_state(state: Mapping[str, object]) -> Mapping[str, object]:
    # PCG64DXSM state is already composed of Python scalars in supported NumPy.
    return json.loads(json.dumps(state, default=lambda value: int(value)))


def atomic_save_npy(path: Path, array: np.ndarray) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
    os.replace(str(temporary), str(path))


def save_placement_batch(path: Path, placement_batches: Sequence[Sequence[PlacedBlock]],
                         arms: Mapping[str, Arm]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw_handle:
        gzip_handle = gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0)
        handle = io.TextIOWrapper(gzip_handle, newline="")
        fields = ["replicate", "block_id", "source_arm", "destination_arm", "chromosome",
                  "block_start0", "block_end0", "components"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for placements in placement_batches:
            for placed in placements:
                writer.writerow({
                    "replicate": placed.replicate,
                    "block_id": placed.block_id,
                    "source_arm": placed.source_arm,
                    "destination_arm": placed.arm,
                    "chromosome": arms[placed.arm].chromosome,
                    "block_start0": placed.start,
                    "block_end0": placed.end,
                    "components": ";".join("%s:%d:%d" % item for item in placed.components),
                })
        handle.flush()
        handle.detach()
        gzip_handle.close()
    os.replace(str(temporary), str(path))


def load_placement_batch(path: Path) -> Iterator[List[PlacedBlock]]:
    current: List[PlacedBlock] = []
    current_replicate: Optional[int] = None
    with gzip.open(str(path), "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            replicate = int(row["replicate"])
            if current_replicate is not None and replicate != current_replicate:
                yield current
                current = []
            current_replicate = replicate
            components = tuple((name, int(offset), int(width))
                               for name, offset, width in (item.split(":") for item in row["components"].split(";")))
            current.append(PlacedBlock(replicate, row["block_id"], row["source_arm"],
                                       row["destination_arm"], int(row["block_start0"]), components))
    if current:
        yield current


def batch_prefix(assignment: str, collection: str, statistic: str, start: int, end: int) -> str:
    return "%s.%s.%s.%09d-%09d" % (assignment, collection, statistic, start, end)


def generate_batches(output: Path, sampler: RegionSampler, genome: GenomeIndex,
                     collections: Sequence[Collection], assignments: Sequence[str],
                     rng: np.random.Generator, start: int, stop: int, batch_size: int,
                     manifest: MutableMapping[str, object]) -> None:
    batches_dir = output / "batches"
    batches_dir.mkdir(exist_ok=True)
    count_dtype = np.uint16 if len(genome.loci) <= np.iinfo(np.uint16).max else np.uint32
    breadth_dtype = np.uint8 if len(genome.arms) <= np.iinfo(np.uint8).max else np.uint16
    cursor = start
    while cursor < stop:
        end = min(stop, cursor + batch_size)
        started = time.perf_counter()
        placements_batch: List[List[PlacedBlock]] = []
        stats_by_assignment: Dict[str, List[AssignmentStats]] = {name: [] for name in assignments}
        for replicate in range(cursor + 1, end + 1):
            placements = sampler.sample(rng, replicate)
            intervals = placed_intervals(placements, sampler.arms)
            placements_batch.append(placements)
            for assignment in assignments:
                stats_by_assignment[assignment].append(
                    calculate_stats(genome, collections, intervals, assignment)
                )
        placement_name = "placements.%09d-%09d.tsv.gz" % (cursor + 1, end)
        save_placement_batch(batches_dir / placement_name, placements_batch, sampler.arms)
        for assignment, values in stats_by_assignment.items():
            atomic_save_npy(batches_dir / (batch_prefix(assignment, "all", "burden", cursor + 1, end) + ".npy"),
                            np.asarray([item.burden for item in values], dtype=np.int32))
            for collection in collections:
                name = collection.name
                atomic_save_npy(batches_dir / (batch_prefix(assignment, name, "counts", cursor + 1, end) + ".npy"),
                                np.stack([item.counts[name] for item in values]).astype(count_dtype))
                atomic_save_npy(batches_dir / (batch_prefix(assignment, name, "denominators", cursor + 1, end) + ".npy"),
                                np.asarray([item.denominators[name] for item in values], dtype=np.int32))
                atomic_save_npy(batches_dir / (batch_prefix(assignment, name, "breadth", cursor + 1, end) + ".npy"),
                                np.stack([item.breadth[name] for item in values]).astype(breadth_dtype))
        elapsed = time.perf_counter() - started
        batch_record = {"start": cursor + 1, "end": end, "placements": placement_name,
                        "elapsed_seconds": elapsed}
        manifest.setdefault("batches", []).append(batch_record)
        manifest["completed_permutations"] = end
        manifest["rng_state"] = jsonable_rng_state(rng.bit_generator.state)
        manifest["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        temporary = output / "run_manifest.json.tmp"
        temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.replace(str(temporary), str(output / "run_manifest.json"))
        cursor = end


def load_batch_array(output: Path, assignment: str, collection: str,
                     statistic: str, record: Mapping[str, object]) -> np.ndarray:
    name = batch_prefix(assignment, collection, statistic,
                        int(record["start"]), int(record["end"])) + ".npy"
    return np.load(str(output / "batches" / name), mmap_mode="r", allow_pickle=False)


def concatenate_batches(output: Path, manifest: Mapping[str, object], assignment: str,
                        collection: str, statistic: str) -> np.ndarray:
    arrays = [np.asarray(load_batch_array(output, assignment, collection, statistic, record))
              for record in manifest["batches"]]
    if not arrays:
        return np.empty(0)
    return np.concatenate(arrays, axis=0)


def concatenate_batch_slice(output: Path, manifest: Mapping[str, object], assignment: str,
                            collection: str, statistic: str, start: int, end: int) -> np.ndarray:
    """Read only a bounded term-column slice from each memory-mapped batch."""
    arrays = [np.asarray(load_batch_array(output, assignment, collection, statistic, record)[:, start:end])
              for record in manifest["batches"]]
    if not arrays:
        return np.empty((0, end - start))
    return np.concatenate(arrays, axis=0)


def summarize_null(values: np.ndarray) -> Tuple[float, float, float, float]:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return math.nan, math.nan, math.nan, math.nan
    q025, median, q975 = np.quantile(values, [0.025, 0.5, 0.975])
    return float(np.mean(values)), float(median), float(q025), float(q975)


def analyze_results(output: Path, manifest: Mapping[str, object], genome: GenomeIndex,
                    blocks: Sequence[Block], arms: Mapping[str, Arm], collections: Sequence[Collection],
                    assignments: Sequence[str]) -> Tuple[List[Mapping[str, object]], List[Mapping[str, object]]]:
    observed_intervals = placed_intervals(observed_placements(blocks), arms)
    observed = {assignment: calculate_stats(genome, collections, observed_intervals, assignment)
                for assignment in assignments}
    permutations = int(manifest["completed_permutations"])
    burden_rows: List[Mapping[str, object]] = []
    result_rows: List[MutableMapping[str, object]] = []
    global_maxima = np.full(permutations, -np.inf, dtype=float)

    for assignment in assignments:
        null_burden = concatenate_batches(output, manifest, assignment, "all", "burden")
        exceed, raw_p, degenerate = empirical_p(observed[assignment].burden, null_burden,
                                                 zero_is_one=True)
        ci_lower, ci_upper = clopper_pearson(exceed, permutations)
        mean, median, q025, q975 = summarize_null(null_burden)
        burden_rows.append({
            "schema_version": SCHEMA_VERSION, "mode": manifest["mode"], "assignment": assignment,
            "run_class": "final" if permutations >= 99_999 else "pilot_nonreportable",
            "reportable_inference": int(permutations >= 99_999),
            "statistic": "all_locus_burden", "observed": observed[assignment].burden,
            "null_mean": mean, "null_median": median, "null_q025": q025, "null_q975": q975,
            "difference_from_median": observed[assignment].burden - median,
            "smoothed_burden_ratio": (observed[assignment].burden + 0.5) / (mean + 0.5),
            "exceedances": exceed, "permutations": permutations, "raw_p": raw_p,
            "mc_ci_lower": ci_lower, "mc_ci_upper": ci_upper,
            "null_variance": float(np.var(null_burden)), "non_informative": int(degenerate),
            "extension_required": int(permutations < 99_999 or exceed < 100),
            "extension_target": (999_999 if (permutations < 99_999 or exceed < 100)
                                 and permutations < 999_999 else ""),
        })
        for collection in collections:
            denominators = concatenate_batches(output, manifest, assignment, collection.name, "denominators")
            null_denominator_median = float(np.median(denominators)) if denominators.size else math.nan
            collection_rows: Dict[str, List[MutableMapping[str, object]]] = {
                statistic: [] for statistic in STATISTICS
            }
            observed_arm_sets: List[Set[str]] = [set() for _ in range(collection.n_terms)]
            for locus_index in observed[assignment].selected:
                locus = genome.loci[int(locus_index)]
                for term_index in collection.locus_terms[int(locus_index)]:
                    observed_arm_sets[int(term_index)].add(locus.arm)
            family_maxima = {
                statistic: np.full(permutations, -np.inf, dtype=float) for statistic in STATISTICS
            }
            # At 100,000 permutations, 64 float64 term columns occupy 51.2 MB.
            # Count/breadth batches remain memory mapped and only these slices
            # are materialized, keeping GO-scale peak memory bounded.
            term_chunk_size = 64
            for term_start in range(0, collection.n_terms, term_chunk_size):
                term_end = min(collection.n_terms, term_start + term_chunk_size)
                counts = concatenate_batch_slice(output, manifest, assignment, collection.name,
                                                 "counts", term_start, term_end)
                breadth = concatenate_batch_slice(output, manifest, assignment, collection.name,
                                                  "breadth", term_start, term_end)
                observed_counts = observed[assignment].counts[collection.name][term_start:term_end]
                observed_breadth = observed[assignment].breadth[collection.name][term_start:term_end]
                observed_denominator = observed[assignment].denominators[collection.name]
                matrices = {
                    "copy_burden": counts.astype(float),
                    "composition": np.divide(
                        counts, denominators[:, None], out=np.zeros_like(counts, dtype=float),
                        where=denominators[:, None] > 0),
                    "breadth": breadth.astype(float),
                }
                observed_values = {
                    "copy_burden": observed_counts.astype(float),
                    "composition": (observed_counts.astype(float) / observed_denominator
                                    if observed_denominator else np.zeros(term_end - term_start)),
                    "breadth": observed_breadth.astype(float),
                }
                for statistic in STATISTICS:
                    matrix = matrices[statistic]
                    obs_values = observed_values[statistic]
                    z_observed, z_null, z_degenerate = standardized_values(obs_values, matrix)
                    if z_null.shape[1]:
                        family_maxima[statistic] = np.maximum(
                            family_maxima[statistic], np.max(z_null, axis=1))
                    for local_index in range(term_end - term_start):
                        term_index = term_start + local_index
                        term_id = collection.term_ids[term_index]
                        null_values = matrix[:, local_index]
                        obs_value = float(obs_values[local_index])
                        exceed, raw_p, degenerate = empirical_p(
                            obs_value, null_values, zero_is_one=(obs_value == 0))
                        ci_lower, ci_upper = clopper_pearson(exceed, permutations)
                        mean, median, q025, q975 = summarize_null(null_values)
                        count_obs = int(observed_counts[local_index])
                        count_null = counts[:, local_index]
                        count_null_median = float(np.median(count_null))
                        cor, effect = composition_effect(
                            count_obs, observed_denominator, count_null_median,
                            null_denominator_median)
                        row: MutableMapping[str, object] = {
                            "schema_version": SCHEMA_VERSION, "mode": manifest["mode"],
                            "run_class": ("final" if permutations >= 99_999 else
                                          "pilot_nonreportable"),
                            "reportable_inference": int(permutations >= 99_999),
                            "assignment": assignment, "collection": collection.name,
                            "statistic": statistic, "term_id": term_id,
                            "term_name": collection.term_names[term_index],
                            "genome_loci": int(collection.genome_locus_counts[term_index]),
                            "genome_arms": int(collection.genome_arm_counts[term_index]),
                            "observed": obs_value, "observed_copy_count": count_obs,
                            "observed_denominator": observed_denominator,
                            "observed_proportion": (count_obs / observed_denominator
                                                    if observed_denominator else 0.0),
                            "observed_arms": ",".join(sorted(
                                observed_arm_sets[term_index], key=canonical_arm_key)),
                            "null_mean": mean, "null_median": median, "null_q025": q025,
                            "null_q975": q975, "null_count_median": count_null_median,
                            "null_denominator_median": null_denominator_median,
                            "smoothed_count_ratio": ((count_obs + 0.5) /
                                                     (float(np.mean(count_null)) + 0.5)),
                            "composition_odds_ratio": cor,
                            "composition_log2_effect": effect,
                            "null_zero_denominator_rate": float(np.mean(denominators == 0)),
                            "exceedances": exceed, "permutations": permutations,
                            "raw_p": raw_p, "mc_ci_lower": ci_lower,
                            "mc_ci_upper": ci_upper,
                            "null_variance": float(np.var(null_values)),
                            "non_informative": int(degenerate or z_degenerate[local_index]),
                            "z_observed": float(z_observed[local_index]),
                        }
                        collection_rows[statistic].append(row)

            for statistic in STATISTICS:
                family_rows = collection_rows[statistic]
                maxima = family_maxima[statistic]
                if assignment == "midpoint" and manifest["mode"] == "primary" and family_rows:
                    global_maxima = np.maximum(global_maxima, maxima)
                p_values = [float(row["raw_p"]) for row in family_rows]
                upper_bounds = [float(row["mc_ci_upper"]) for row in family_rows]
                bh = bh_adjust(p_values)
                bh_upper = bh_adjust(upper_bounds)
                by = by_adjust(p_values)
                for index, row in enumerate(family_rows):
                    row["bh_q"] = float(bh[index])
                    row["bh_upper_ci_q"] = float(bh_upper[index])
                    row["by_q"] = float(by[index])
                    row["family_hypotheses"] = len(family_rows)
                    if row["non_informative"] or float(row["observed"]) == 0:
                        max_exceed = permutations
                    else:
                        max_exceed = int(np.count_nonzero(maxima >= float(row["z_observed"])))
                    max_ci = clopper_pearson(max_exceed, permutations)
                    row["collection_maxT_exceedances"] = max_exceed
                    row["collection_maxT_p"] = (max_exceed + 1.0) / (permutations + 1.0)
                    row["collection_maxT_ci_lower"] = max_ci[0]
                    row["collection_maxT_ci_upper"] = max_ci[1]
                    null_always_zero = (float(row["observed"]) > 0 and
                                        float(row["null_mean"]) == 0)
                    near_threshold = any(
                        0.04 <= float(value) <= 0.06
                        for value in (row["raw_p"], row["bh_q"],
                                      row["collection_maxT_p"]))
                    row["extension_required"] = int(
                        permutations < 99_999 or int(row["exceedances"]) < 100 or
                        int(row["collection_maxT_exceedances"]) < 100 or
                        null_always_zero or near_threshold)
                    row["extension_target"] = (999_999 if row["extension_required"] and
                                                permutations < 999_999 else "")
                result_rows.extend(family_rows)

    # Global primary maxT contains midpoint term hypotheses only.  Sensitivity
    # assignments and non-primary backgrounds remain explicitly excluded.
    if "midpoint" in assignments and manifest["mode"] == "primary":
        primary_rows = [row for row in result_rows if row["assignment"] == "midpoint"]
        for row in result_rows:
            if row["assignment"] != "midpoint":
                row["global_maxT_exceedances"] = ""
                row["global_maxT_p"] = ""
                row["global_maxT_ci_lower"] = ""
                row["global_maxT_ci_upper"] = ""
                continue
            if row["non_informative"] or float(row["observed"]) == 0:
                exceed = permutations
            else:
                exceed = int(np.count_nonzero(global_maxima >= float(row["z_observed"])))
            row["global_maxT_exceedances"] = exceed
            row["global_maxT_p"] = (exceed + 1.0) / (permutations + 1.0)
            bounds = clopper_pearson(exceed, permutations)
            row["global_maxT_ci_lower"], row["global_maxT_ci_upper"] = bounds
            if exceed < 100 or 0.04 <= float(row["global_maxT_p"]) <= 0.06:
                row["extension_required"] = 1
                row["extension_target"] = 999_999 if permutations < 999_999 else ""
    else:
        for row in result_rows:
            row["global_maxT_exceedances"] = ""
            row["global_maxT_p"] = ""
            row["global_maxT_ci_lower"] = ""
            row["global_maxT_ci_upper"] = ""
    # The any-overlap assignment is a paired sensitivity using exactly the same
    # placements.  Boundary sensitivity is defined from composition effect and
    # attached to every midpoint statistic for the term.
    overlap_effect = {
        (str(row["collection"]), str(row["term_id"])): float(row["composition_log2_effect"])
        for row in result_rows
        if row["assignment"] == "overlap" and row["statistic"] == "composition"
    }
    midpoint_effect = {
        (str(row["collection"]), str(row["term_id"])): float(row["composition_log2_effect"])
        for row in result_rows
        if row["assignment"] == "midpoint" and row["statistic"] == "composition"
    }
    for row in result_rows:
        row["boundary_sensitive"] = ""
        if row["assignment"] != "midpoint":
            continue
        key = (str(row["collection"]), str(row["term_id"]))
        if key not in overlap_effect or key not in midpoint_effect:
            continue
        primary, sensitivity = midpoint_effect[key], overlap_effect[key]
        direction_change = ((primary > 0 and sensitivity <= 0) or
                            (primary < 0 and sensitivity >= 0))
        effect_change = abs(sensitivity - primary) > 0.5 * max(abs(primary), 0.25)
        row["boundary_sensitive"] = int(direction_change or effect_change)
    return burden_rows, result_rows


def group_concentration(groups: Sequence[str]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    if not groups:
        return None, None, None
    counts = Counter(groups)
    total = sum(counts.values())
    fractions = [count / total for count in counts.values()]
    herfindahl = sum(value * value for value in fractions)
    return max(fractions), herfindahl, 1.0 / herfindahl


def driver_rows(result_rows: Sequence[Mapping[str, object]], observed_stats: AssignmentStats,
                collection_by_name: Mapping[str, Collection], genome: GenomeIndex,
                family_groups: np.ndarray, identity_groups: np.ndarray) -> Tuple[List[Mapping[str, object]], List[Mapping[str, object]]]:
    diagnostic_keys: Set[Tuple[str, str]] = set()
    for row in result_rows:
        if row["assignment"] != "midpoint":
            continue
        global_value = row.get("global_maxT_p", "")
        if float(row["bh_q"]) <= 0.10 or (global_value != "" and float(global_value) <= 0.10):
            diagnostic_keys.add((str(row["collection"]), str(row["term_id"])))
    summaries: List[Mapping[str, object]] = []
    groups_out: List[Mapping[str, object]] = []
    for collection_name, term_id in sorted(diagnostic_keys):
        collection = collection_by_name[collection_name]
        term_index = collection.term_ids.index(term_id)
        selected = [int(index) for index in observed_stats.selected
                    if term_index in collection.locus_terms[int(index)]]
        family_values = [str(family_groups[index]) for index in selected]
        identity_values = [str(identity_groups[index]) for index in selected]
        family_fraction, family_hhi, family_neff = group_concentration(family_values)
        identity_fraction, identity_hhi, identity_neff = group_concentration(identity_values)
        summaries.append({
            "collection": collection_name, "term_id": term_id, "observed_copies": len(selected),
            "distinct_locus_ids": len({genome.loci[index].locus_id for index in selected}),
            "distinct_symbols": len({genome.loci[index].gene_name for index in selected}),
            "distinct_arms": len({genome.loci[index].arm for index in selected}),
            "largest_family_fraction": family_fraction, "family_herfindahl": family_hhi,
            "family_effective_number": family_neff,
            "largest_identity_fraction": identity_fraction, "identity_herfindahl": identity_hhi,
            "identity_effective_number": identity_neff,
            "multiple_annotations_same_locus": int(any(
                len(collection.locus_terms[index]) > 1 for index in selected)),
        })
        for grouping, values in (("symbol", [genome.loci[index].gene_name for index in selected]),
                                 ("arm", [genome.loci[index].arm for index in selected]),
                                 ("family", family_values), ("identity", identity_values),
                                 ("locus_id", [genome.loci[index].locus_id for index in selected])):
            for group, count in sorted(Counter(values).items(), key=lambda item: (-item[1], item[0])):
                groups_out.append({"collection": collection_name, "term_id": term_id,
                                   "grouping": grouping, "group_id": group, "copy_count": count})
    return summaries, groups_out


def _term_metrics(genome: GenomeIndex, collection: Collection, selected: np.ndarray,
                  term_index: int, grouping: np.ndarray, excluded_group: str) -> Tuple[int, int, float, int]:
    retained = [int(index) for index in selected if str(grouping[int(index)]) != excluded_group]
    denominator = sum(1 for index in retained if collection.annotated[index])
    term_loci = [index for index in retained if term_index in collection.locus_terms[index]]
    count = len(term_loci)
    breadth = len({genome.loci[index].arm for index in term_loci})
    proportion = count / denominator if denominator else 0.0
    return count, denominator, proportion, breadth


def leave_one_diagnostics(output: Path, manifest: Mapping[str, object],
                          result_rows: Sequence[Mapping[str, object]], genome: GenomeIndex,
                          arms: Mapping[str, Arm], blocks: Sequence[Block],
                          collections: Mapping[str, Collection], family_groups: np.ndarray,
                          identity_groups: np.ndarray) -> Tuple[List[Mapping[str, object]], Dict[Tuple[str, str], Mapping[str, object]]]:
    """Paired leave-one-family/identity influence analyses.

    Every removal is applied to the same loci in the observation and all saved
    randomized region sets.  No placement is resampled.
    """
    diagnostic_keys: Set[Tuple[str, str]] = set()
    base_effects: Dict[Tuple[str, str], float] = {}
    for row in result_rows:
        if row["assignment"] != "midpoint":
            continue
        key = (str(row["collection"]), str(row["term_id"]))
        if row["statistic"] == "composition":
            base_effects[key] = float(row["composition_log2_effect"])
        global_value = row.get("global_maxT_p", "")
        if float(row["bh_q"]) <= 0.10 or (global_value != "" and float(global_value) <= 0.10):
            diagnostic_keys.add(key)
    if not diagnostic_keys:
        return [], {}

    observed_intervals = placed_intervals(observed_placements(blocks), arms)
    observed_selected, _ = genome.select(observed_intervals, "midpoint")
    requests: List[MutableMapping[str, object]] = []
    for key in sorted(diagnostic_keys):
        collection_name, term_id = key
        collection = collections[collection_name]
        term_index = collection.term_ids.index(term_id)
        contributing = [int(index) for index in observed_selected
                        if term_index in collection.locus_terms[int(index)]]
        for grouping_name, grouping in (("family", family_groups), ("identity", identity_groups)):
            counts = Counter(str(grouping[index]) for index in contributing)
            ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
            if len(ordered) <= 20:
                chosen = ordered
            else:
                chosen = []
                covered = 0
                for item in ordered:
                    chosen.append(item)
                    covered += item[1]
                    if covered >= 0.8 * len(contributing):
                        break
            for rank, (group_id, observed_contribution) in enumerate(chosen, 1):
                count, denominator, proportion, breadth = _term_metrics(
                    genome, collection, observed_selected, term_index, grouping, group_id)
                request: MutableMapping[str, object] = {
                    "collection": collection_name, "term_id": term_id,
                    "grouping": grouping_name, "group_id": group_id,
                    "contribution_rank": rank, "observed_removed_copies": observed_contribution,
                    "observed_copy_count": count, "observed_denominator": denominator,
                    "observed_composition": proportion, "observed_breadth": breadth,
                    "_collection": collection, "_term_index": term_index, "_grouping": grouping,
                    "_null_counts": [], "_null_denominators": [], "_null_composition": [],
                    "_null_breadth": [],
                }
                requests.append(request)

    # Scan each auditable placement exactly once and update every requested
    # removal.  This is intentionally deferred until diagnostic terms are known.
    observed_replicates = 0
    for batch in manifest["batches"]:
        placement_path = output / "batches" / str(batch["placements"])
        for placements in load_placement_batch(placement_path):
            selected, _ = genome.select(placed_intervals(placements, arms), "midpoint")
            observed_replicates += 1
            for request in requests:
                count, denominator, proportion, breadth = _term_metrics(
                    genome, request["_collection"], selected, int(request["_term_index"]),
                    request["_grouping"], str(request["group_id"]))
                request["_null_counts"].append(count)
                request["_null_denominators"].append(denominator)
                request["_null_composition"].append(proportion)
                request["_null_breadth"].append(breadth)
    if observed_replicates != int(manifest["completed_permutations"]):
        raise RuntimeError("saved placement count disagrees with run manifest during leave-one analysis")

    rows: List[Mapping[str, object]] = []
    flags: Dict[Tuple[str, str], Mapping[str, object]] = {}
    for request in requests:
        null_counts = np.asarray(request.pop("_null_counts"), dtype=float)
        null_denominators = np.asarray(request.pop("_null_denominators"), dtype=float)
        null_composition = np.asarray(request.pop("_null_composition"), dtype=float)
        null_breadth = np.asarray(request.pop("_null_breadth"), dtype=float)
        request.pop("_collection")
        request.pop("_term_index")
        request.pop("_grouping")
        count = int(request["observed_copy_count"])
        denominator = int(request["observed_denominator"])
        breadth = int(request["observed_breadth"])
        count_exceed, count_p, _ = empirical_p(count, null_counts, zero_is_one=(count == 0))
        composition_exceed, composition_p, _ = empirical_p(
            float(request["observed_composition"]), null_composition,
            zero_is_one=(float(request["observed_composition"]) == 0))
        breadth_exceed, breadth_p, _ = empirical_p(breadth, null_breadth, zero_is_one=(breadth == 0))
        cor, effect = composition_effect(count, denominator, float(np.median(null_counts)),
                                         float(np.median(null_denominators)))
        key = (str(request["collection"]), str(request["term_id"]))
        base_effect = base_effects.get(key, 0.0)
        direction_reversed = ((base_effect > 0 and effect <= 0) or
                              (base_effect < 0 and effect >= 0))
        effect_reduced_half = (abs(base_effect) > 0 and abs(effect) <= 0.5 * abs(base_effect))
        p_above = max(count_p, composition_p, breadth_p) > 0.05
        sensitive = direction_reversed or breadth <= 1 or p_above or effect_reduced_half
        request.update({
            "copy_burden_exceedances": count_exceed, "copy_burden_raw_p": count_p,
            "composition_exceedances": composition_exceed, "composition_raw_p": composition_p,
            "breadth_exceedances": breadth_exceed, "breadth_raw_p": breadth_p,
            "null_count_median": float(np.median(null_counts)),
            "null_denominator_median": float(np.median(null_denominators)),
            "null_composition_median": float(np.median(null_composition)),
            "null_breadth_median": float(np.median(null_breadth)),
            "composition_odds_ratio": cor, "composition_log2_effect": effect,
            "direction_reversed": int(direction_reversed),
            "effect_reduced_at_least_half": int(effect_reduced_half),
            "raw_p_above_0.05": int(p_above), "single_driver_sensitive": int(sensitive),
            "permutations": observed_replicates,
        })
        rows.append(request)
        existing = flags.get(key, {"single_driver_sensitive": 0,
                                   "largest_family_breadth": "",
                                   "largest_identity_breadth": ""})
        existing = dict(existing)
        existing["single_driver_sensitive"] = max(int(existing["single_driver_sensitive"]), int(sensitive))
        if int(request["contribution_rank"]) == 1:
            existing["largest_%s_breadth" % request["grouping"]] = breadth
        flags[key] = existing
    return rows, flags


def parse_named_path(text: str) -> Tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("expected NAME=PATH")
    name, path = text.split("=", 1)
    if not name or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise argparse.ArgumentTypeError("collection name must match [A-Za-z0-9_.-]+")
    return name, Path(path)


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def immutable_configuration(args: argparse.Namespace, files: Mapping[str, Path],
                            collections: Sequence[Tuple[str, Path]]) -> Mapping[str, object]:
    checksums = {name: {"path": str(path.resolve()), "bytes": path.stat().st_size,
                        "sha256": sha256_file(path)} for name, path in files.items() if path is not None}
    return {
        "schema_version": SCHEMA_VERSION, "engine_version": ENGINE_VERSION,
        "mode": args.mode, "assignments": args.assignment,
        "seed": args.seed, "bit_generator": "PCG64DXSM", "spawn_key": [0],
        "threads": 1,
        "batch_size": args.batch_size, "min_candidates": args.min_candidates,
        "min_term_loci": args.min_term_loci, "min_term_arms": args.min_term_arms,
        "mask_enumeration_limit": args.mask_enumeration_limit,
        "mask_policy": "frozen_bed" if args.mask is not None else "empty",
        "pilot_allowed": bool(args.allow_pilot),
        "driver_diagnostics_deferred": bool(args.skip_driver_diagnostics),
        "collections": [{"name": name, "path": str(path.resolve()), "sha256": sha256_file(path),
                         "bytes": path.stat().st_size} for name, path in collections],
        "inputs": checksums,
    }


def run(args: argparse.Namespace) -> None:
    started = time.perf_counter()
    output = args.output
    named_collections = [parse_named_path(value) for value in args.terms]
    if len({name for name, _path in named_collections}) != len(named_collections):
        raise ValueError("collection names must be unique within a run")
    for _name, path in named_collections:
        if not path.is_file():
            raise FileNotFoundError(path)
    files = {"arms": args.arms, "intervals": args.intervals, "loci": args.loci,
             "mask": args.mask, "family_map": args.family_map,
             "identity_map": args.identity_map, "engine": Path(__file__),
             "statistical_spec": Path(__file__).with_name("STATISTICAL_SPEC.md")}
    for name, path in files.items():
        if path is not None and not path.is_file():
            raise FileNotFoundError("%s: %s" % (name, path))
    config = immutable_configuration(args, files, named_collections)
    manifest_path = output / "run_manifest.json"
    if output.exists() and not args.resume:
        raise FileExistsError("output exists; choose a new directory or pass --resume: %s" % output)
    if args.resume:
        if not manifest_path.is_file():
            raise FileNotFoundError("--resume requires run_manifest.json")
        manifest: MutableMapping[str, object] = json.loads(manifest_path.read_text())
        if manifest["immutable_configuration"] != config:
            raise ValueError("resume configuration differs from the frozen original run")
        if int(manifest["completed_permutations"]) > args.permutations:
            raise ValueError("cannot resume to fewer permutations than already completed")
    else:
        output.mkdir(parents=True)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "immutable_configuration": config,
            "mode": args.mode,
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "completed_permutations": 0,
            "batches": [],
            "command": " ".join(sys.argv),
            "git_commit": git_commit(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "threads": 1,
        }

    arms = load_arms(args.arms)
    intervals = load_intervals(args.intervals, arms)
    loci = load_loci(args.loci, arms)
    genome = GenomeIndex(loci, arms)
    blocks = build_blocks(intervals, arms)
    masks = load_masks(args.mask, arms)
    sampler = RegionSampler(blocks, arms, masks, args.mode, args.min_candidates,
                            args.mask_enumeration_limit)
    collections = [load_collection(name, path, genome, args.min_term_loci, args.min_term_arms)
                   for name, path in named_collections]
    family_groups = load_group_map(args.family_map, genome, "family")
    identity_groups = load_group_map(args.identity_map, genome, "identity")
    assignments = args.assignment

    if not args.resume:
        manifest["target_interval_count"] = len(intervals)
        manifest["target_union_bp"] = sum(end - start for arm in arms for start, end in
                                          coalesce_pairs((row.start, row.end) for row in intervals if row.arm == arm))
        manifest["physical_locus_count"] = len(loci)
        manifest["placement_block_count"] = len(blocks)
        seed_sequence = np.random.SeedSequence(args.seed)
        child = seed_sequence.spawn(1)[0]
        rng = np.random.Generator(np.random.PCG64DXSM(child))
        manifest["rng_state"] = jsonable_rng_state(rng.bit_generator.state)
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        rng = np.random.Generator(np.random.PCG64DXSM())
        rng.bit_generator.state = manifest["rng_state"]

    candidate_fields = ["mode", "block_id", "source_arm", "destination_arm", "candidate_count",
                        "range_count", "explicit", "observed_start_is_candidate"]
    write_tsv(output / "candidate_spaces.tsv", candidate_fields, sampler.candidate_rows())
    block_rows = []
    for block in blocks:
        block_rows.append({"block_id": block.block_id, "source_arm": block.source_arm,
                           "source_start0": block.source_start, "source_end0": block.source_end,
                           "span_bp": block.span, "stratum": block.stratum,
                           "component_count": len(block.components),
                           "components": ";".join("%s:%d:%d" % item for item in block.components)})
    write_tsv(output / "placement_blocks.tsv",
              ["block_id", "source_arm", "source_start0", "source_end0", "span_bp", "stratum",
               "component_count", "components"], block_rows)
    filtered = [row for collection in collections for row in collection.filtered_rows]
    write_tsv(output / "term_filtering.tsv",
              ["collection", "term_id", "term_name", "genome_loci", "genome_arms", "status"], filtered)

    completed = int(manifest["completed_permutations"])
    if args.permutations > completed:
        generate_batches(output, sampler, genome, collections, assignments, rng, completed,
                         args.permutations, args.batch_size, manifest)
        manifest = json.loads(manifest_path.read_text())

    burden, results = analyze_results(output, manifest, genome, blocks, arms, collections, assignments)
    burden_fields = list(burden[0].keys()) if burden else []
    write_tsv(output / "burden_results.tsv", burden_fields, burden)
    result_fields = list(results[0].keys()) if results else [
        "schema_version", "mode", "assignment", "collection", "statistic", "term_id", "term_name"
    ]
    write_tsv(output / "term_results.tsv", result_fields, results)

    observed_midpoint = calculate_stats(genome, collections,
                                        placed_intervals(observed_placements(blocks), arms), "midpoint")
    if args.skip_driver_diagnostics:
        # Large final runs may finalize the prespecified driver analysis in a
        # separate, exact cached-subtraction pass. This switch does not alter
        # placement generation, cached null arrays, or term inference.
        summaries, group_counts, leave_rows, leave_flags = [], [], [], {}
    else:
        summaries, group_counts = driver_rows(results, observed_midpoint,
                                              {item.name: item for item in collections}, genome,
                                              family_groups, identity_groups)
        leave_rows, leave_flags = leave_one_diagnostics(
            output, manifest, results, genome, arms, blocks,
            {item.name: item for item in collections}, family_groups, identity_groups)
    for summary in summaries:
        key = (str(summary["collection"]), str(summary["term_id"]))
        values = leave_flags.get(key, {})
        summary["breadth_after_largest_family"] = values.get("largest_family_breadth", "")
        summary["breadth_after_largest_identity"] = values.get("largest_identity_breadth", "")
        summary["single_driver_sensitive"] = values.get("single_driver_sensitive", "")
    summary_fields = list(summaries[0].keys()) if summaries else [
        "collection", "term_id", "observed_copies", "distinct_locus_ids", "distinct_symbols",
        "distinct_arms", "largest_family_fraction", "family_herfindahl", "family_effective_number",
        "largest_identity_fraction", "identity_herfindahl", "identity_effective_number",
        "multiple_annotations_same_locus", "breadth_after_largest_family",
        "breadth_after_largest_identity", "single_driver_sensitive"]
    write_tsv(output / "driver_summary.tsv", summary_fields, summaries)
    group_fields = list(group_counts[0].keys()) if group_counts else [
        "collection", "term_id", "grouping", "group_id", "copy_count"]
    write_tsv(output / "driver_group_counts.tsv", group_fields, group_counts)
    leave_fields = list(leave_rows[0].keys()) if leave_rows else [
        "collection", "term_id", "grouping", "group_id", "contribution_rank",
        "observed_removed_copies", "observed_copy_count", "observed_denominator",
        "observed_composition", "observed_breadth", "copy_burden_exceedances",
        "copy_burden_raw_p", "composition_exceedances", "composition_raw_p",
        "breadth_exceedances", "breadth_raw_p", "null_count_median",
        "null_denominator_median", "null_composition_median", "null_breadth_median",
        "composition_odds_ratio", "composition_log2_effect", "direction_reversed",
        "effect_reduced_at_least_half", "raw_p_above_0.05", "single_driver_sensitive",
        "permutations"]
    write_tsv(output / "leave_one_sensitivity.tsv", leave_fields, leave_rows)

    if args.mode == "adjacent":
        comparator = sampler.immediate_adjacent()
        comparator_rows = []
        comparator_intervals = placed_intervals(comparator, arms)
        for assignment in assignments:
            values = calculate_stats(genome, collections, comparator_intervals, assignment)
            comparator_rows.append({"assignment": assignment, "all_locus_burden": values.burden,
                                    "block_count": len(comparator), "p_value": "not_applicable"})
        write_tsv(output / "immediate_adjacent_comparator.tsv",
                  ["assignment", "all_locus_burden", "block_count", "p_value"], comparator_rows)

    prior_qc = manifest.get("placement_qc_counts", {})
    all_qc_arms = set(prior_qc) | set(sampler.accepted_by_arm) | set(sampler.rejected_by_arm)
    cumulative_qc = {}
    qc_rows = []
    for arm_name in sorted(all_qc_arms, key=canonical_arm_key):
        previous = prior_qc.get(arm_name, {})
        accepted = int(previous.get("accepted", 0)) + sampler.accepted_by_arm[arm_name]
        rejected = int(previous.get("rejected", 0)) + sampler.rejected_by_arm[arm_name]
        cumulative_qc[arm_name] = {"accepted": accepted, "rejected": rejected}
        qc_rows.append({"arm": arm_name, "accepted_joint_draws": accepted,
                        "rejected_joint_draws": rejected,
                        "acceptance_rate": accepted / (accepted + rejected) if accepted + rejected else math.nan})
    write_tsv(output / "placement_qc.tsv",
              ["arm", "accepted_joint_draws", "rejected_joint_draws", "acceptance_rate"], qc_rows)
    manifest["placement_qc_counts"] = cumulative_qc
    manifest["elapsed_seconds_total_last_invocation"] = time.perf_counter() - started
    manifest["output_files"] = sorted(path.name for path in output.iterdir() if path.is_file())
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def build_parser() -> argparse.ArgumentParser:
    directory = Path(__file__).with_name("analysis_ready")
    parser = argparse.ArgumentParser(
        description="Copy-number-aware CHM13 PHR coherent-region permutation engine")
    parser.add_argument("--arms", type=Path, default=directory / "chm13_arm_summary.tsv")
    parser.add_argument("--intervals", type=Path, default=directory / "chm13_phr_intervals.tsv")
    parser.add_argument("--loci", type=Path, default=directory / "chm13_gene_loci.tsv.gz")
    parser.add_argument("--terms", action="append", default=[], metavar="NAME=PATH",
                        help="frozen locus-to-term TSV; repeat for arbitrary declared collections")
    parser.add_argument("--family-map", type=Path)
    parser.add_argument("--identity-map", type=Path)
    parser.add_argument("--mask", type=Path, help="optional frozen exclusion-mask BED")
    parser.add_argument("--mode", choices=sorted(DEFAULT_SEEDS), default="primary")
    parser.add_argument("--assignment", choices=("midpoint", "overlap"), action="append",
                        help="repeat to calculate both; default: midpoint and overlap")
    parser.add_argument("--seed", type=int, help="explicit decimal master seed")
    parser.add_argument("--permutations", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--min-candidates", type=int, default=100)
    parser.add_argument("--min-term-loci", type=int, default=5)
    parser.add_argument("--min-term-arms", type=int, default=2)
    parser.add_argument("--mask-enumeration-limit", type=int, default=10_000_000)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--skip-driver-diagnostics", action="store_true",
                        help="write empty driver tables for exact separate finalization")
    parser.add_argument("--allow-pilot", action="store_true",
                        help="permit <99,999 permutations; outputs are marked non-reportable")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.permutations < 1 or args.batch_size < 1:
        parser.error("--permutations and --batch-size must be positive")
    if args.permutations < 99_999 and not args.allow_pilot:
        parser.error("fewer than 99,999 permutations requires --allow-pilot; pilot p-values are non-reportable")
    if args.min_candidates < 1 or args.min_term_loci < 1 or args.min_term_arms < 1:
        parser.error("minimum candidate/term thresholds must be positive")
    if args.assignment is None:
        args.assignment = ["midpoint", "overlap"]
    elif len(set(args.assignment)) != len(args.assignment):
        parser.error("duplicate --assignment")
    if args.seed is None:
        parser.error("--seed is required (prespecified %s seed is %d)" %
                     (args.mode, DEFAULT_SEEDS[args.mode]))
    try:
        run(args)
    except (ValueError, FileNotFoundError, FileExistsError, RuntimeError) as exc:
        parser.exit(2, "COPY_engine.py: error: %s\n" % exc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
