#!/usr/bin/env python3
"""
Summarize completed native odgi untangle BEDPE/PAF merge-distance outputs.

Run from the repository root:

    python3 scripts/pedigree/untangle_native_merge_tracts.py

Large native odgi intermediates are generated outside git by the Slurm runner:

    scripts/pedigree/run_untangle_native_merge_tracts.sbatch

The committed output is a compact summary table in scripts/pedigree/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
ODGI = Path("/home/erikg/.guix-profile/bin/odgi")
GRAPH = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/pggb/"
    "washu.1Mb.telo_500kb_trimmed.fa.gz.6e0e250.11fba48.13f423a.smooth.final.og"
)
UNTANGLE_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle")
EXTERNAL_OUT = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_slurm")
COMMUNITY_TSV = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/flanking/similarity/hprcv2.1Mb.flanking.communities.tsv"
)
PREVIOUS_SUMMARY = REPO_ROOT / "scripts" / "pedigree" / "untangle_multimap_tract_summary.tsv"

OUT_SUMMARY = REPO_ROOT / "scripts" / "pedigree" / "untangle_native_merge_summary.tsv"
REPORT = REPO_ROOT / "paper_prep" / "_brainstorming" / "pedigree_native_untangle_merge.md"
SWEEPGA = Path("/home/erikg/.cargo/bin/sweepga")
SWEEPGA_REPO = Path("/moosefs/erikg/sweepga")
SWEEPGA_COMMIT = "018e4ce49d2c125820e0ac50dc5feaa02d423683"
SWEEPGA_TEST_DIR = Path("/moosefs/erikg/phrs/pedigree_native_untangle_agent2556_sweepga_test")

CONVERSION_RANGE = (22, 95)
CROSSOVER_RANGE = (318, 688)
NEAR_CROSSOVER_RANGE = (159, 1376)

PAIRINGS = {
    "PAN027_vs_PAN010": {
        "child": "PAN027",
        "parent": "PAN010",
        "child_hap": "1",
        "label": "PAN027 maternal (hap1) vs PAN010 (mother)",
    },
    "PAN027_vs_PAN011": {
        "child": "PAN027",
        "parent": "PAN011",
        "child_hap": "2",
        "label": "PAN027 paternal (hap2) vs PAN011 (father)",
    },
    "PAN028_vs_PAN027": {
        "child": "PAN028",
        "parent": "PAN027",
        "child_hap": "1",
        "label": "PAN028 maternal (hap1) vs PAN027 (mother)",
    },
}


@dataclass(frozen=True)
class NativeFile:
    pair: str
    merge_dist: int
    n_best_emit: int
    fmt: str
    path: Path

    @property
    def dataset(self) -> str:
        return f"native_{self.fmt}.e50000.m{self.merge_dist}.n{self.n_best_emit}.j0.8"


@dataclass(frozen=True)
class Hit:
    sample: str
    hap: str
    chrom: str
    chrarm: str
    community: str
    score: float
    nth_best: int
    strand: str
    target_name: str
    target_start: int
    target_end: int

    @property
    def hap_key(self) -> str:
        return f"{self.sample}#{self.hap}:{self.chrarm}"


@dataclass
class Segment:
    dataset: str
    source_format: str
    setting: str
    pair: str
    label: str
    query_name: str
    query_chr: str
    query_chrarm: str
    start: int
    end: int
    hits: list[Hit]

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def hap_keys(self) -> set[str]:
        return {hit.hap_key for hit in self.hits}

    @property
    def arm_keys(self) -> set[str]:
        return {hit.chrarm for hit in self.hits}

    @property
    def communities(self) -> set[str]:
        return {hit.community for hit in self.hits if hit.community != "unknown"}


@dataclass
class Tract:
    dataset: str
    source_format: str
    setting: str
    pair: str
    label: str
    query_name: str
    query_chr: str
    query_chrarm: str
    start: int
    end: int
    segments: list[Segment]

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--external-out", type=Path, default=EXTERNAL_OUT)
    parser.add_argument("--emit-n-best", type=int, default=4, help="native odgi -n value to emit")
    parser.add_argument("--analysis-n-best", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument(
        "--merge-dist",
        type=int,
        nargs="+",
        default=[1000],
        help="completed native odgi merge distances to summarize; pass 0 1000 for a full audit",
    )
    parser.add_argument("--min-jaccard", type=float, default=0.8)
    return parser.parse_args()


def parse_path_info(path_name: str) -> tuple[str, str, str, str] | None:
    match = re.search(r"^(PAN\d+)#([12])#.*_(chr(?:[0-9]+|X|Y))_([pq]arm)$", path_name)
    if not match:
        return None
    sample, hap, chrom, arm = match.groups()
    return sample, hap, chrom, chrom + ("p" if arm == "parm" else "q")


def open_text(path: Path):
    if path.suffix == ".gz":
        with path.open("rb") as test:
            if test.read(2) != b"\x1f\x8b":
                return path.open(errors="replace")
        return gzip.open(path, "rt")
    return path.open(errors="replace")


def load_communities(path: Path) -> dict[str, str]:
    communities: dict[str, str] = {}
    if not path.exists():
        return communities
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            arm = (row.get("arm") or "").replace("_", "")
            if arm.startswith("chr"):
                communities[arm] = row.get("community", "unknown") or "unknown"
    return communities


def native_files(args: argparse.Namespace) -> list[NativeFile]:
    files: list[NativeFile] = []
    for pair in PAIRINGS:
        for merge_dist in args.merge_dist:
            for fmt in ("bedpe", "paf"):
                suffix = "bedpe.gz" if fmt == "bedpe" else "paf.gz"
                path = (
                    args.external_out
                    / f"{pair}.e50000.m{merge_dist}.j{args.min_jaccard}.n{args.emit_n_best}.{suffix}"
                )
                files.append(NativeFile(pair, merge_dist, args.emit_n_best, fmt, path))
    return files


def odgi_command(native: NativeFile, args: argparse.Namespace, query_list: Path) -> list[str]:
    target_list = UNTANGLE_DIR / f"targets_{native.pair}.txt"
    cmd = [
        str(ODGI),
        "untangle",
        "-t",
        "${SLURM_CPUS_PER_TASK:-24}",
        "-i",
        str(GRAPH),
        "-Q",
        str(query_list),
        "-R",
        str(target_list),
        "-e",
        "50000",
        "-m",
        str(native.merge_dist),
        "-j",
        str(args.min_jaccard),
        "-n",
        str(native.n_best_emit),
    ]
    if native.fmt == "paf":
        cmd.append("-p")
    return cmd


def paf_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_native_file(native: NativeFile, communities: dict[str, str]) -> list[Segment]:
    grouped: dict[tuple[str, int, int], list[Hit]] = defaultdict(list)
    with open_text(native.path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if native.fmt == "bedpe":
                if len(fields) < 10:
                    continue
                query_name, q_start, q_end = fields[0], int(fields[1]), int(fields[2])
                target_name, t_start, t_end = fields[3], int(fields[4]), int(fields[5])
                score = float(fields[6])
                strand = fields[7]
                nth_best = int(float(fields[9]))
            else:
                if len(fields) < 12:
                    continue
                tags = paf_tags(fields)
                query_name, q_start, q_end = fields[0], int(fields[2]), int(fields[3])
                target_name, t_start, t_end = fields[5], int(fields[7]), int(fields[8])
                score = float(tags.get("jc", tags.get("sc", "nan")))
                strand = fields[4]
                nth_best = int(tags.get("nb", "1"))
            query_info = parse_path_info(query_name)
            target_info = parse_path_info(target_name)
            if query_info is None or target_info is None:
                continue
            sample, hap, chrom, chrarm = target_info
            grouped[(query_name, q_start, q_end)].append(
                Hit(
                    sample=sample,
                    hap=hap,
                    chrom=chrom,
                    chrarm=chrarm,
                    community=communities.get(chrarm, "unknown"),
                    score=score,
                    nth_best=nth_best,
                    strand=strand,
                    target_name=target_name,
                    target_start=t_start,
                    target_end=t_end,
                )
            )
    segments: list[Segment] = []
    info = PAIRINGS[native.pair]
    for (query_name, q_start, q_end), hits in grouped.items():
        query_info = parse_path_info(query_name)
        if query_info is None:
            continue
        _, _, query_chr, query_chrarm = query_info
        segments.append(
            Segment(
                dataset=native.dataset,
                source_format=native.fmt,
                setting="",
                pair=native.pair,
                label=info["label"],
                query_name=query_name,
                query_chr=query_chr,
                query_chrarm=query_chrarm,
                start=q_start,
                end=q_end,
                hits=sorted(hits, key=lambda hit: (hit.nth_best, -hit.score, hit.sample, hit.hap, hit.chrarm)),
            )
        )
    segments.sort(key=lambda segment: (segment.query_name, segment.start, segment.end))
    return segments


def iter_native_segments(native: NativeFile, communities: dict[str, str]):
    """Yield native segments from odgi output without retaining the full file.

    odgi untangle writes rows sorted by query path and query interval. The
    streaming parser relies on that order and flushes one query interval at a
    time, which keeps m0 PAF/BEDPE parsing bounded in memory.
    """
    current_key: tuple[str, int, int] | None = None
    current_hits: list[Hit] = []

    def flush():
        if current_key is None or not current_hits:
            return None
        query_name, q_start, q_end = current_key
        query_info = parse_path_info(query_name)
        if query_info is None:
            return None
        _, _, query_chr, query_chrarm = query_info
        info = PAIRINGS[native.pair]
        return Segment(
            dataset=native.dataset,
            source_format=native.fmt,
            setting="",
            pair=native.pair,
            label=info["label"],
            query_name=query_name,
            query_chr=query_chr,
            query_chrarm=query_chrarm,
            start=q_start,
            end=q_end,
            hits=sorted(current_hits, key=lambda hit: (hit.nth_best, -hit.score, hit.sample, hit.hap, hit.chrarm)),
        )

    with open_text(native.path) as handle:
        for line in handle:
            parsed = parse_native_line(native, line, communities)
            if parsed is None:
                continue
            key, hit = parsed
            if current_key is not None and key != current_key:
                segment = flush()
                if segment is not None:
                    yield segment
                current_hits = []
            current_key = key
            current_hits.append(hit)
    segment = flush()
    if segment is not None:
        yield segment


def parse_native_line(
    native: NativeFile, line: str, communities: dict[str, str]
) -> tuple[tuple[str, int, int], Hit] | None:
    if not line.strip() or line.startswith("#"):
        return None
    fields = line.rstrip("\n").split("\t")
    try:
        if native.fmt == "bedpe":
            if len(fields) < 10:
                return None
            query_name, q_start, q_end = fields[0], int(fields[1]), int(fields[2])
            target_name, t_start, t_end = fields[3], int(fields[4]), int(fields[5])
            score = float(fields[6])
            strand = fields[7]
            nth_best = int(float(fields[9]))
        else:
            if len(fields) < 12:
                return None
            tags = paf_tags(fields)
            query_name, q_start, q_end = fields[0], int(fields[2]), int(fields[3])
            target_name, t_start, t_end = fields[5], int(fields[7]), int(fields[8])
            score = float(tags.get("jc", tags.get("sc", "nan")))
            strand = fields[4]
            nth_best = int(tags.get("nb", "1"))
    except ValueError:
        return None
    if parse_path_info(query_name) is None:
        return None
    target_info = parse_path_info(target_name)
    if target_info is None:
        return None
    sample, hap, chrom, chrarm = target_info
    return (
        (query_name, q_start, q_end),
        Hit(
            sample=sample,
            hap=hap,
            chrom=chrom,
            chrarm=chrarm,
            community=communities.get(chrarm, "unknown"),
            score=score,
            nth_best=nth_best,
            strand=strand,
            target_name=target_name,
            target_start=t_start,
            target_end=t_end,
        ),
    )


def classify_hits(hits: Iterable[Hit]) -> str:
    hits = list(hits)
    if not hits:
        return "unresolved_no_call"
    hap_keys = {hit.hap_key for hit in hits}
    arms = {hit.chrarm for hit in hits}
    communities = {hit.community for hit in hits if hit.community != "unknown"}
    if len(hap_keys) == 1:
        return "unique_donor_haplotype"
    if len(arms) == 1:
        return "unique_donor_arm_multi_haplotype"
    if len(communities) == 1:
        return "same_community_ambiguous"
    return "cross_community_ambiguous"


def compatible(left: Segment, right: Segment) -> bool:
    if left.hap_keys & right.hap_keys:
        return True
    if left.arm_keys & right.arm_keys:
        return True
    if left.communities & right.communities:
        return True
    return False


def project_segments(
    native: NativeFile,
    segments: list[Segment],
    analysis_n: int,
    min_jaccard: float,
) -> list[Segment]:
    setting = f"{native.fmt}_m{native.merge_dist}_n{analysis_n}"
    projected: list[Segment] = []
    for segment in segments:
        hits = [hit for hit in segment.hits if hit.nth_best <= analysis_n and hit.score >= min_jaccard]
        if not hits:
            continue
        projected.append(
            Segment(
                dataset=native.dataset,
                source_format=native.fmt,
                setting=setting,
                pair=segment.pair,
                label=segment.label,
                query_name=segment.query_name,
                query_chr=segment.query_chr,
                query_chrarm=segment.query_chrarm,
                start=segment.start,
                end=segment.end,
                hits=hits,
            )
        )
    return projected


def project_one_segment(
    native: NativeFile,
    segment: Segment,
    analysis_n: int,
    min_jaccard: float,
) -> Segment | None:
    hits = [hit for hit in segment.hits if hit.nth_best <= analysis_n and hit.score >= min_jaccard]
    if not hits:
        return None
    return Segment(
        dataset=native.dataset,
        source_format=native.fmt,
        setting=f"{native.fmt}_m{native.merge_dist}_n{analysis_n}",
        pair=segment.pair,
        label=segment.label,
        query_name=segment.query_name,
        query_chr=segment.query_chr,
        query_chrarm=segment.query_chrarm,
        start=segment.start,
        end=segment.end,
        hits=hits,
    )


def merge_native_segments(segments: list[Segment]) -> list[Tract]:
    tracts: list[Tract] = []
    by_query: dict[str, list[Segment]] = defaultdict(list)
    for segment in segments:
        by_query[segment.query_name].append(segment)
    for query_name in sorted(by_query):
        current: Tract | None = None
        for segment in sorted(by_query[query_name], key=lambda item: (item.start, item.end)):
            if current is None:
                current = start_tract(segment)
                continue
            if segment.start <= current.end and compatible(current.segments[-1], segment):
                current.end = max(current.end, segment.end)
                current.segments.append(segment)
            else:
                tracts.append(current)
                current = start_tract(segment)
        if current is not None:
            tracts.append(current)
    return tracts


def merge_native_segments_stream(segments: Iterable[Segment]) -> list[Tract]:
    tracts: list[Tract] = []
    current: Tract | None = None
    current_query: str | None = None
    def keep(tract: Tract) -> None:
        if interchromosomal(tract):
            tracts.append(tract)

    for segment in segments:
        if current is None or current_query != segment.query_name:
            if current is not None:
                keep(current)
            current = start_tract(segment)
            current_query = segment.query_name
            continue
        if segment.start <= current.end and compatible(current.segments[-1], segment):
            current.end = max(current.end, segment.end)
            current.segments.append(segment)
        else:
            keep(current)
            current = start_tract(segment)
            current_query = segment.query_name
    if current is not None:
        keep(current)
    return tracts


def start_tract(segment: Segment) -> Tract:
    return Tract(
        dataset=segment.dataset,
        source_format=segment.source_format,
        setting=segment.setting,
        pair=segment.pair,
        label=segment.label,
        query_name=segment.query_name,
        query_chr=segment.query_chr,
        query_chrarm=segment.query_chrarm,
        start=segment.start,
        end=segment.end,
        segments=[segment],
    )


def interchromosomal(tract: Tract) -> bool:
    return any(hit.chrom != tract.query_chr for segment in tract.segments for hit in segment.hits)


def tract_resolvability(tract: Tract) -> str:
    return classify_hits(hit for segment in tract.segments for hit in segment.hits)


def donor_meta(tract: Tract) -> dict[str, str]:
    hits = [hit for segment in tract.segments for hit in segment.hits]
    scores = [hit.score for hit in hits]
    strands = sorted({hit.strand for hit in hits})
    return {
        "donor_haplotypes": ",".join(sorted({hit.hap_key for hit in hits})),
        "donor_arms": ",".join(sorted({hit.chrarm for hit in hits})),
        "donor_communities": ",".join(sorted({hit.community for hit in hits})),
        "strands": ",".join(strands),
        "min_score": fmt_float(min(scores) if scores else None),
        "max_score": fmt_float(max(scores) if scores else None),
    }


def fmt_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def quantile(values: list[int], prob: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return float(ordered[0])
    h = (len(ordered) - 1) * prob
    lo = int(h)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] * (1 - (h - lo)) + ordered[hi] * (h - lo)


def in_range(value: int, bounds: tuple[int, int]) -> bool:
    return bounds[0] <= value <= bounds[1]


def proportion(n: int, d: int) -> str:
    return f"{n / d:.6f}" if d else "NA"


def write_segments(segments: list[Segment], path: Path) -> None:
    rows = []
    for segment in segments:
        for hit in segment.hits:
            rows.append(
                {
                    "dataset": segment.dataset,
                    "setting": segment.setting,
                    "source_format": segment.source_format,
                    "pair": segment.pair,
                    "query_name": segment.query_name,
                    "query_chrarm": segment.query_chrarm,
                    "query_start": segment.start,
                    "query_end": segment.end,
                    "target_name": hit.target_name,
                    "target_start": hit.target_start,
                    "target_end": hit.target_end,
                    "strand": hit.strand,
                    "nth_best": hit.nth_best,
                    "score": fmt_float(hit.score),
                    "donor_haplotype": hit.hap_key,
                    "donor_arm": hit.chrarm,
                    "donor_community": hit.community,
                    "segment_resolvability": classify_hits(segment.hits),
                }
            )
    write_tsv(rows, path)


def write_tracts(tracts: list[Tract], path: Path) -> None:
    rows = []
    for tract in tracts:
        meta = donor_meta(tract)
        row = {
            "dataset": tract.dataset,
            "setting": tract.setting,
            "source_format": tract.source_format,
            "pair": tract.pair,
            "label": tract.label,
            "query_name": tract.query_name,
            "query_chrarm": tract.query_chrarm,
            "start": tract.start,
            "end": tract.end,
            "length_bp": tract.length,
            "n_segments": len(tract.segments),
            "resolvability": tract_resolvability(tract),
            "is_interchromosomal": interchromosomal(tract),
        }
        row.update(meta)
        rows.append(row)
    write_tsv(rows, path)


def summary_rows(tracts: list[Tract], previous: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    by_setting: dict[tuple[str, str], list[Tract]] = defaultdict(list)
    for tract in tracts:
        by_setting[(tract.source_format, tract.setting)].append(tract)

    n_by_key = {
        (tract.source_format, tract.setting): len([t for t in values if interchromosomal(t)])
        for (tract.source_format, tract.setting), values in by_setting.items()
    }
    rows = []
    for (source_format, setting), values in sorted(by_setting.items()):
        inter = [tract for tract in values if interchromosomal(tract)]
        lengths = [tract.length for tract in inter]
        n = len(inter)
        classes = Counter(tract_resolvability(tract) for tract in inter)
        merge_dist = int(re.search(r"_m(\d+)_", setting).group(1))
        analysis_n = int(re.search(r"_n(\d+)$", setting).group(1))
        counterpart = f"{source_format}_m{0 if merge_dist == 1000 else 1000}_n{analysis_n}"
        m0_n = n_by_key.get((source_format, f"{source_format}_m0_n{analysis_n}"))
        m1000_n = n_by_key.get((source_format, f"{source_format}_m1000_n{analysis_n}"))
        native_merged_vs_m0 = (
            str(max(0, (m0_n or 0) - n)) if merge_dist == 1000 and m0_n is not None else "NA"
        )
        native_split_vs_m0 = (
            str(max(0, n - (m0_n or 0))) if merge_dist == 1000 and m0_n is not None else "NA"
        )
        previous_default = previous.get("m1000_default", {})
        rows.append(
            {
                "source_format": source_format,
                "setting": setting,
                "n_tracts": str(n),
                "native_m0_tracts_for_same_n": str(m0_n) if m0_n is not None else "NA",
                "native_m1000_tracts_for_same_n": str(m1000_n) if m1000_n is not None else "NA",
                "native_m1000_merged_vs_m0_n": native_merged_vs_m0,
                "native_m1000_split_vs_m0_n": native_split_vs_m0,
                "counterpart_setting": counterpart,
                "prior_posthoc_m1000_default_n": previous_default.get("n_tracts", "NA"),
                "prior_posthoc_m1000_bridge1kb_n": previous.get("m1000_bridge1kb", {}).get("n_tracts", "NA"),
                "unique_donor_haplotype_n": str(classes["unique_donor_haplotype"]),
                "unique_donor_arm_multi_haplotype_n": str(classes["unique_donor_arm_multi_haplotype"]),
                "same_community_ambiguous_n": str(classes["same_community_ambiguous"]),
                "cross_community_ambiguous_n": str(classes["cross_community_ambiguous"]),
                "unresolved_no_call_n": str(classes["unresolved_no_call"]),
                "ambiguous_tract_n": str(
                    classes["unique_donor_arm_multi_haplotype"]
                    + classes["same_community_ambiguous"]
                    + classes["cross_community_ambiguous"]
                    + classes["unresolved_no_call"]
                ),
                "min_bp": fmt_float(float(min(lengths)) if lengths else None),
                "q25_bp": fmt_float(quantile(lengths, 0.25)),
                "median_bp": fmt_float(quantile(lengths, 0.50)),
                "q75_bp": fmt_float(quantile(lengths, 0.75)),
                "max_bp": fmt_float(float(max(lengths)) if lengths else None),
                "conversion_22_95_n": str(sum(in_range(length, CONVERSION_RANGE) for length in lengths)),
                "conversion_22_95_prop": proportion(sum(in_range(length, CONVERSION_RANGE) for length in lengths), n),
                "crossover_318_688_n": str(sum(in_range(length, CROSSOVER_RANGE) for length in lengths)),
                "crossover_318_688_prop": proportion(sum(in_range(length, CROSSOVER_RANGE) for length in lengths), n),
                "near_crossover_159_1376_n": str(sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths)),
                "near_crossover_159_1376_prop": proportion(
                    sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths), n
                ),
            }
        )
    return rows


def summarize_setting(
    source_format: str,
    setting: str,
    tracts: list[Tract],
    previous: dict[str, dict[str, str]],
) -> dict[str, str]:
    lengths = [tract.length for tract in tracts]
    n = len(tracts)
    classes = Counter(tract_resolvability(tract) for tract in tracts)
    previous_default = previous.get("m1000_default", {})
    return {
        "source_format": source_format,
        "setting": setting,
        "n_tracts": str(n),
        "native_m0_tracts_for_same_n": "NA",
        "native_m1000_tracts_for_same_n": "NA",
        "native_m1000_merged_vs_m0_n": "NA",
        "native_m1000_split_vs_m0_n": "NA",
        "counterpart_setting": "NA",
        "prior_posthoc_m1000_default_n": previous_default.get("n_tracts", "NA"),
        "prior_posthoc_m1000_bridge1kb_n": previous.get("m1000_bridge1kb", {}).get("n_tracts", "NA"),
        "unique_donor_haplotype_n": str(classes["unique_donor_haplotype"]),
        "unique_donor_arm_multi_haplotype_n": str(classes["unique_donor_arm_multi_haplotype"]),
        "same_community_ambiguous_n": str(classes["same_community_ambiguous"]),
        "cross_community_ambiguous_n": str(classes["cross_community_ambiguous"]),
        "unresolved_no_call_n": str(classes["unresolved_no_call"]),
        "ambiguous_tract_n": str(
            classes["unique_donor_arm_multi_haplotype"]
            + classes["same_community_ambiguous"]
            + classes["cross_community_ambiguous"]
            + classes["unresolved_no_call"]
        ),
        "min_bp": fmt_float(float(min(lengths)) if lengths else None),
        "q25_bp": fmt_float(quantile(lengths, 0.25)),
        "median_bp": fmt_float(quantile(lengths, 0.50)),
        "q75_bp": fmt_float(quantile(lengths, 0.75)),
        "max_bp": fmt_float(float(max(lengths)) if lengths else None),
        "conversion_22_95_n": str(sum(in_range(length, CONVERSION_RANGE) for length in lengths)),
        "conversion_22_95_prop": proportion(sum(in_range(length, CONVERSION_RANGE) for length in lengths), n),
        "crossover_318_688_n": str(sum(in_range(length, CROSSOVER_RANGE) for length in lengths)),
        "crossover_318_688_prop": proportion(sum(in_range(length, CROSSOVER_RANGE) for length in lengths), n),
        "near_crossover_159_1376_n": str(sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths)),
        "near_crossover_159_1376_prop": proportion(
            sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths), n
        ),
    }


def empty_summary_agg() -> dict[str, object]:
    return {"n": 0, "lengths": [], "classes": Counter()}


def update_summary_agg(agg: dict[str, object], tracts: list[Tract]) -> None:
    lengths = agg["lengths"]
    classes = agg["classes"]
    assert isinstance(lengths, list)
    assert isinstance(classes, Counter)
    agg["n"] = int(agg["n"]) + len(tracts)
    lengths.extend(tract.length for tract in tracts)
    classes.update(tract_resolvability(tract) for tract in tracts)


def update_summary_agg_one(agg: dict[str, object], tract: Tract) -> None:
    lengths = agg["lengths"]
    classes = agg["classes"]
    assert isinstance(lengths, list)
    assert isinstance(classes, Counter)
    agg["n"] = int(agg["n"]) + 1
    lengths.append(tract.length)
    classes.update([tract_resolvability(tract)])


def consume_native_stream(
    native: NativeFile,
    args: argparse.Namespace,
    communities: dict[str, str],
    summary_aggs: dict[tuple[str, str], dict[str, object]],
) -> None:
    current: dict[int, Tract | None] = {analysis_n: None for analysis_n in args.analysis_n_best}

    def keep(analysis_n: int, tract: Tract | None) -> None:
        if tract is None or not interchromosomal(tract):
            return
        setting = f"{native.fmt}_m{native.merge_dist}_n{analysis_n}"
        update_summary_agg_one(summary_aggs[(native.fmt, setting)], tract)

    for segment in iter_native_segments(native, communities):
        for analysis_n in args.analysis_n_best:
            projected = project_one_segment(native, segment, analysis_n, args.min_jaccard)
            if projected is None:
                continue
            tract = current[analysis_n]
            if tract is None or tract.query_name != projected.query_name:
                keep(analysis_n, tract)
                current[analysis_n] = start_tract(projected)
            elif projected.start <= tract.end and compatible(tract.segments[-1], projected):
                tract.end = max(tract.end, projected.end)
                tract.segments.append(projected)
            else:
                keep(analysis_n, tract)
                current[analysis_n] = start_tract(projected)

    for analysis_n, tract in current.items():
        keep(analysis_n, tract)


def summarize_agg(
    source_format: str,
    setting: str,
    agg: dict[str, object],
    previous: dict[str, dict[str, str]],
) -> dict[str, str]:
    lengths = agg["lengths"]
    classes = agg["classes"]
    assert isinstance(lengths, list)
    assert isinstance(classes, Counter)
    n = int(agg["n"])
    previous_default = previous.get("m1000_default", {})
    return {
        "source_format": source_format,
        "setting": setting,
        "n_tracts": str(n),
        "native_m0_tracts_for_same_n": "NA",
        "native_m1000_tracts_for_same_n": "NA",
        "native_m1000_merged_vs_m0_n": "NA",
        "native_m1000_split_vs_m0_n": "NA",
        "counterpart_setting": "NA",
        "prior_posthoc_m1000_default_n": previous_default.get("n_tracts", "NA"),
        "prior_posthoc_m1000_bridge1kb_n": previous.get("m1000_bridge1kb", {}).get("n_tracts", "NA"),
        "unique_donor_haplotype_n": str(classes["unique_donor_haplotype"]),
        "unique_donor_arm_multi_haplotype_n": str(classes["unique_donor_arm_multi_haplotype"]),
        "same_community_ambiguous_n": str(classes["same_community_ambiguous"]),
        "cross_community_ambiguous_n": str(classes["cross_community_ambiguous"]),
        "unresolved_no_call_n": str(classes["unresolved_no_call"]),
        "ambiguous_tract_n": str(
            classes["unique_donor_arm_multi_haplotype"]
            + classes["same_community_ambiguous"]
            + classes["cross_community_ambiguous"]
            + classes["unresolved_no_call"]
        ),
        "min_bp": fmt_float(float(min(lengths)) if lengths else None),
        "q25_bp": fmt_float(quantile(lengths, 0.25)),
        "median_bp": fmt_float(quantile(lengths, 0.50)),
        "q75_bp": fmt_float(quantile(lengths, 0.75)),
        "max_bp": fmt_float(float(max(lengths)) if lengths else None),
        "conversion_22_95_n": str(sum(in_range(length, CONVERSION_RANGE) for length in lengths)),
        "conversion_22_95_prop": proportion(sum(in_range(length, CONVERSION_RANGE) for length in lengths), n),
        "crossover_318_688_n": str(sum(in_range(length, CROSSOVER_RANGE) for length in lengths)),
        "crossover_318_688_prop": proportion(sum(in_range(length, CROSSOVER_RANGE) for length in lengths), n),
        "near_crossover_159_1376_n": str(sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths)),
        "near_crossover_159_1376_prop": proportion(
            sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths), n
        ),
    }


def add_native_count_comparisons(rows: list[dict[str, str]]) -> None:
    by_key = {(row["source_format"], row["setting"]): row for row in rows}
    counts = {(row["source_format"], row["setting"]): int(row["n_tracts"]) for row in rows}
    for row in rows:
        source_format = row["source_format"]
        setting = row["setting"]
        merge_dist = int(re.search(r"_m(\d+)_", setting).group(1))
        analysis_n = int(re.search(r"_n(\d+)$", setting).group(1))
        m0_setting = f"{source_format}_m0_n{analysis_n}"
        m1000_setting = f"{source_format}_m1000_n{analysis_n}"
        m0_n = counts.get((source_format, m0_setting))
        m1000_n = counts.get((source_format, m1000_setting))
        row["native_m0_tracts_for_same_n"] = str(m0_n) if m0_n is not None else "NA"
        row["native_m1000_tracts_for_same_n"] = str(m1000_n) if m1000_n is not None else "NA"
        row["counterpart_setting"] = m1000_setting if merge_dist == 0 else m0_setting
        if merge_dist == 1000 and m0_n is not None:
            n = int(row["n_tracts"])
            row["native_m1000_merged_vs_m0_n"] = str(max(0, m0_n - n))
            row["native_m1000_split_vs_m0_n"] = str(max(0, n - m0_n))


def load_previous_summary(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open() as handle:
        return {row["setting"]: row for row in csv.DictReader(handle, delimiter="\t")}


def write_tsv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_visual(tracts: list[Tract]) -> None:
    VIS_DIR.mkdir(parents=True, exist_ok=True)
    candidates = [
        tract
        for tract in tracts
        if tract.source_format == "paf"
        and tract.setting == "paf_m1000_n4"
        and interchromosomal(tract)
        and tract.length <= 150000
    ]
    by_class: dict[str, list[Tract]] = defaultdict(list)
    for tract in candidates:
        by_class[tract_resolvability(tract)].append(tract)
    selected: list[Tract] = []
    for cls in [
        "unique_donor_haplotype",
        "unique_donor_arm_multi_haplotype",
        "same_community_ambiguous",
        "cross_community_ambiguous",
    ]:
        selected.extend(sorted(by_class.get(cls, []), key=lambda tract: (-len(tract.segments), -tract.length))[:4])
    selected = selected[:16]
    rows = []
    for tract in selected:
        meta = donor_meta(tract)
        for segment in tract.segments:
            rows.append(
                {
                    "setting": tract.setting,
                    "pair": tract.pair,
                    "query_name": tract.query_name,
                    "tract_start": tract.start,
                    "tract_end": tract.end,
                    "segment_start": segment.start,
                    "segment_end": segment.end,
                    "resolvability": tract_resolvability(tract),
                    "segment_resolvability": classify_hits(segment.hits),
                    "donor_arms": meta["donor_arms"],
                    "donor_communities": meta["donor_communities"],
                    "min_score": meta["min_score"],
                    "max_score": meta["max_score"],
                }
            )
    write_tsv(rows, VIS_TSV)
    write_visual_svg(rows, VIS_SVG)


def write_visual_svg(rows: list[dict[str, object]], path: Path) -> None:
    if not rows:
        path.write_text("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"80\"></svg>\n")
        return
    colors = {
        "unique_donor_haplotype": "#1b9e77",
        "unique_donor_arm_multi_haplotype": "#7570b3",
        "same_community_ambiguous": "#e6ab02",
        "cross_community_ambiguous": "#d95f02",
        "unresolved_no_call": "#666666",
    }
    keys = []
    for row in rows:
        key = (row["setting"], row["pair"], row["query_name"], row["tract_start"], row["tract_end"])
        if key not in keys:
            keys.append(key)
    width = 1200
    left = 95
    right = 470
    row_h = 28
    top = 52
    height = top + row_h * len(keys) + 86
    min_x = min(int(row["tract_start"]) for row in rows)
    max_x = max(int(row["tract_end"]) for row in rows)
    span = max(1, max_x - min_x)

    def sx(value: object) -> float:
        return left + (int(value) - min_x) / span * (width - left - right)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,sans-serif;font-size:11px}.title{font-size:16px;font-weight:bold}.label{font-size:9px}</style>",
        f'<text class="title" x="{left}" y="24">Native odgi PAF m1000 n4 representative multimapping tracts</text>',
    ]
    for y, key in enumerate(keys):
        group = [
            row
            for row in rows
            if (row["setting"], row["pair"], row["query_name"], row["tract_start"], row["tract_end"]) == key
        ]
        cy = top + y * row_h
        parts.append(
            f'<line x1="{sx(key[3]):.1f}" x2="{sx(key[4]):.1f}" y1="{cy}" y2="{cy}" '
            'stroke="#222" stroke-width="1" opacity="0.45"/>'
        )
        for row in group:
            cls = str(row["segment_resolvability"])
            x = sx(row["segment_start"])
            w = max(1, sx(row["segment_end"]) - x)
            parts.append(
                f'<rect x="{x:.1f}" y="{cy-6}" width="{w:.1f}" height="12" '
                f'fill="{colors.get(cls, "#999999")}"><title>{html.escape(cls)}</title></rect>'
            )
        label = f"{key[1]} {str(key[2]).split('#')[-1]} {int(key[4]) - int(key[3])} bp"
        parts.append(f'<text class="label" x="{width-right+12}" y="{cy+4}">{html.escape(label)}</text>')
    legend_y = height - 62
    legend_x = left
    for cls, color in colors.items():
        parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="12" height="12" fill="{color}"/>')
        parts.append(f'<text class="label" x="{legend_x+17}" y="{legend_y+10}">{html.escape(cls)}</text>')
        legend_x += 225
        if legend_x > width - 250:
            legend_x = left
            legend_y += 18
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def command_manifest(args: argparse.Namespace) -> list[str]:
    lines: list[str] = []
    for pair in PAIRINGS:
        query_list = args.external_out / f"queries_{pair}.txt"
        for native in native_files(args):
            if native.pair != pair:
                continue
            cmd = odgi_command(native, args, query_list)
            lines.append(" ".join(cmd) + f" | gzip -c > {native.path}")
    return lines


def count_lines(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open(errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def base_summary_row(source_format: str, setting: str, n_tracts: str) -> dict[str, str]:
    return {
        "source_format": source_format,
        "setting": setting,
        "n_tracts": n_tracts,
        "native_m0_tracts_for_same_n": "NA",
        "native_m1000_tracts_for_same_n": "NA",
        "native_m1000_merged_vs_m0_n": "NA",
        "native_m1000_split_vs_m0_n": "NA",
        "counterpart_setting": "NA",
        "prior_posthoc_m1000_default_n": "NA",
        "prior_posthoc_m1000_bridge1kb_n": "NA",
        "unique_donor_haplotype_n": "NA",
        "unique_donor_arm_multi_haplotype_n": "NA",
        "same_community_ambiguous_n": "NA",
        "cross_community_ambiguous_n": "NA",
        "unresolved_no_call_n": "NA",
        "ambiguous_tract_n": "NA",
        "min_bp": "NA",
        "q25_bp": "NA",
        "median_bp": "NA",
        "q75_bp": "NA",
        "max_bp": "NA",
        "conversion_22_95_n": "NA",
        "conversion_22_95_prop": "NA",
        "crossover_318_688_n": "NA",
        "crossover_318_688_prop": "NA",
        "near_crossover_159_1376_n": "NA",
        "near_crossover_159_1376_prop": "NA",
    }


def add_sweepga_rows(summary: list[dict[str, str]]) -> None:
    for num_mappings in ("1_many", "2_many", "4_many"):
        path = (
            SWEEPGA_TEST_DIR
            / f"PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.{num_mappings}.paf"
        )
        n = count_lines(path)
        summary.append(
            base_summary_row(
                "sweepga_paf",
                f"PAN028_vs_PAN027_m1000_n4_num_mappings_{num_mappings}",
                str(n) if n is not None else "MISSING",
            )
        )


def sweepga_version() -> str:
    proc = subprocess.run([str(SWEEPGA), "--version"], text=True, capture_output=True, check=False)
    return (proc.stdout or proc.stderr).strip() or "unavailable"


def sweepga_repo_commit() -> str:
    proc = subprocess.run(
        ["git", "-C", str(SWEEPGA_REPO), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
        check=False,
    )
    return proc.stdout.strip() or "unavailable"


def sweepga_status(num_mappings: str) -> tuple[str, str]:
    safe = num_mappings.replace(":", "_")
    stderr = SWEEPGA_TEST_DIR / f"uncompressed_{safe}.stderr"
    out = SWEEPGA_TEST_DIR / f"PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.{safe}.paf"
    gz_stderr = SWEEPGA_TEST_DIR / f"{safe}.stderr"
    status = "accepted" if out.exists() and out.stat().st_size > 0 else "not_available"
    detail = ""
    if stderr.exists():
        text = stderr.read_text(errors="replace").strip().splitlines()
        detail = text[-1] if text else ""
    if gz_stderr.exists():
        gz_text = gz_stderr.read_text(errors="replace")
        if "invalid BGZF header" in gz_text:
            detail = (detail + "; " if detail else "") + "gzip path rejected: invalid BGZF header"
    return status, detail


def write_report(args: argparse.Namespace, commands: list[str], summary: list[dict[str, str]]) -> None:
    def row(setting: str, source_format: str = "paf") -> dict[str, str]:
        for item in summary:
            if item["setting"] == setting and item["source_format"] == source_format:
                return item
        return {}

    paf_m0_n4 = row("paf_m0_n4")
    paf_m1000_n4 = row("paf_m1000_n4")
    bed_m1000_n4 = row("bedpe_m1000_n4", "bedpe")
    version = sweepga_version()
    commit = sweepga_repo_commit()
    sweepga_lines = []
    for num_mappings in ("1:many", "2:many", "4:many"):
        safe = num_mappings.replace(":", "_")
        output = SWEEPGA_TEST_DIR / f"PAN028_vs_PAN027.e50000.m1000.j0.8.n4.uncompressed.{safe}.paf"
        n = count_lines(output)
        status, detail = sweepga_status(num_mappings)
        sweepga_lines.append(
            f"- `{SWEEPGA} --num-mappings {num_mappings} --scaffold-jump 0 "
            f"--output-file {output} {SWEEPGA_TEST_DIR / 'PAN028_vs_PAN027.e50000.m1000.j0.8.n4.native.paf'}`: "
            f"{status}, {n if n is not None else 'NA'} PAF rows"
            + (f" ({detail})" if detail else "")
        )

    text = f"""# Pedigree Native odgi BEDPE/PAF Merge-Distance Decision Record

Primary runnable: `scripts/pedigree/run_untangle_native_merge_tracts.sbatch`.
Parser: `scripts/pedigree/untangle_native_merge_tracts.py`.

Large native odgi intermediates are outside git under:

- `{args.external_out}`

## Were native BEDPE/PAF outputs generated on Slurm?

Yes. Slurm job `1703959` generated the valid rerun outputs in
`{args.external_out}` using 24 CPUs, 96G, and the `workers` partition. The first
worker also ran a direct head-node pass before the Slurm-only constraint was
issued; that earlier pass is not treated as the valid provenance. The committed
sbatch script is now the natural runnable and directly runs both BEDPE and PAF
forms of `odgi untangle`, including `-p` for PAF.

The Slurm grid was three WashU child-parent comparisons, merge distances `0`
and `1000`, `-e 50000`, `-j {args.min_jaccard}`, and `-n {args.emit_n_best}`.
The default compact summary parses the `m1000` files because that is sufficient
for the BEDPE/PAF/sweepGA comparison; pass `--merge-dist 0 1000` for a full
native merge-distance audit. The parser projects odgi's emitted `nb`/`nth.best`
fields to analysis top-N values 1, 2, and 4 without rerunning heavy work.

Command manifest:

```bash
{chr(10).join(commands)}
```

## Does sweepGA accept/filter odgi-emitted PAF directly?

Yes for an uncompressed odgi-emitted PAF. The installed binary is
`{SWEEPGA}` reporting `{version}`, from `/moosefs/erikg/sweepga` commit
`{commit}`. The required commit for this task is `{SWEEPGA_COMMIT}`.

Minimal representative test used
`PAN028_vs_PAN027.e50000.m1000.j0.8.n4.paf.gz` from the Slurm output. Passing
the `.paf.gz` path directly failed with `invalid BGZF header`; after plain
decompression to `.paf`, sweepGA accepted and filtered the native odgi PAF:

{chr(10).join(sweepga_lines)}

The compact comparison is in `scripts/pedigree/untangle_native_merge_summary.tsv`.

## Does native merge or sweepGA filtering justify a later manuscript edit?

No. Native BEDPE/PAF output and sweepGA filtering establish a cleaner
provenance path, but they do not clearly improve tract calls enough to justify
a later manuscript edit from this task alone. For orientation only, native PAF
`m1000/n4` has {paf_m1000_n4.get('n_tracts', 'NA')} parser-counted
interchromosomal tracts and native BEDPE `m1000/n4` has
{bed_m1000_n4.get('n_tracts', 'NA')}. These are methods/provenance counts, not
a conversion-vs-crossover mechanism claim.
"""
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text)


def validate_odgi_help() -> None:
    proc = subprocess.run([str(ODGI), "untangle", "--help"], text=True, capture_output=True, check=True)
    required = ["--paf-output", "--merge-dist", "--n-best", "--min-jaccard", "--cut-every", "BEDPE"]
    missing = [item for item in required if item not in proc.stdout]
    if missing:
        raise RuntimeError(f"odgi untangle --help is missing expected text: {', '.join(missing)}")


def main() -> int:
    args = parse_args()
    validate_odgi_help()
    commands = command_manifest(args)

    missing = [native.path for native in native_files(args) if not native.path.exists()]
    if missing:
        for path in missing:
            print(f"missing native output: {path}", file=sys.stderr)
        print("submit scripts/pedigree/run_untangle_native_merge_tracts.sbatch to generate outputs", file=sys.stderr)
        return 2

    communities = load_communities(COMMUNITY_TSV)
    summary_aggs: dict[tuple[str, str], dict[str, object]] = defaultdict(empty_summary_agg)
    previous = load_previous_summary(PREVIOUS_SUMMARY)
    for native in native_files(args):
        consume_native_stream(native, args, communities, summary_aggs)

    summary = [
        summarize_agg(source_format, setting, agg, previous)
        for (source_format, setting), agg in sorted(summary_aggs.items())
    ]
    add_native_count_comparisons(summary)
    add_sweepga_rows(summary)
    write_tsv(summary, OUT_SUMMARY)
    write_report(args, commands, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
