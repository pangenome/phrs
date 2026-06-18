#!/usr/bin/env python3
"""
Call multimap-aware candidate tracts from WashU odgi untangle BED files.

Run from the repository root:

    python3 scripts/pedigree/untangle_multimap_tracts.py

The caller groups all untangle rows with the same child/query interval into an
equivalence class, retains top-N rows that are tied or near-tied to the best
score, and merges adjacent intervals when their donor equivalence classes are
compatible under a chosen bridge mode.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[2]
UNTANGLE_DIR = Path("/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu/untangle")
COMMUNITY_TSV = Path(
    "/moosefs/guarracino/HPRCv2/PHR_III/flanking/similarity/hprcv2.1Mb.flanking.communities.tsv"
)
PATCH_SUMMARY_TSV = REPO_ROOT / "scripts" / "pedigree" / "patch_tract_length_summary.tsv"
LOWER_SUMMARY_TSV = REPO_ROOT / "scripts" / "pedigree" / "patch_tract_lower_merge_summary.tsv"
OUT_TRACTS = REPO_ROOT / "scripts" / "pedigree" / "untangle_multimap_tracts.tsv"
OUT_SUMMARY = REPO_ROOT / "scripts" / "pedigree" / "untangle_multimap_tract_summary.tsv"
OUT_DIR = REPO_ROOT / "paper_prep" / "_brainstorming" / "pedigree_multimap_tracts"
OUT_MD = REPO_ROOT / "paper_prep" / "_brainstorming" / "pedigree_multimap_tracts.md"

CONVERSION_RANGE = (22, 95)
CROSSOVER_RANGE = (318, 688)
NEAR_CONVERSION_RANGE = (11, 190)
NEAR_CROSSOVER_RANGE = (159, 1376)

PAIRINGS = {
    "PAN027_vs_PAN010": ("PAN027#1", "PAN027 maternal (hap1) vs PAN010 (mother)"),
    "PAN027_vs_PAN011": ("PAN027#2", "PAN027 paternal (hap2) vs PAN011 (father)"),
    "PAN028_vs_PAN027": ("PAN028#1", "PAN028 maternal (hap1) vs PAN027 (mother)"),
}

DEFAULT_RUNS = (
    "m1000_default:1000:4:0.001:0.8:0:arm",
    "m1000_bridge1kb:1000:8:0.002:0.8:1000:community",
)


@dataclass(frozen=True)
class Hit:
    sample: str
    hap: str
    chrom: str
    chrarm: str
    community: str
    score: float
    nth_best: int
    ref_name: str

    @property
    def hap_key(self) -> str:
        return f"{self.sample}#{self.hap}:{self.chrarm}"

    @property
    def arm_key(self) -> str:
        return self.chrarm


@dataclass
class Segment:
    dataset: str
    setting: str
    pair: str
    label: str
    query_name: str
    query_chr: str
    query_chrarm: str
    start: int
    end: int
    best_score: float
    hits: list[Hit]

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def hap_keys(self) -> set[str]:
        return {hit.hap_key for hit in self.hits}

    @property
    def arm_keys(self) -> set[str]:
        return {hit.arm_key for hit in self.hits}

    @property
    def communities(self) -> set[str]:
        return {hit.community for hit in self.hits if hit.community != "unknown"}


@dataclass
class Tract:
    dataset: str
    setting: str
    pair: str
    label: str
    query_name: str
    query_chr: str
    query_chrarm: str
    start: int
    end: int
    segments: list[Segment]
    bridge_bp: int
    bridge_reasons: Counter[str]

    @property
    def length(self) -> int:
        return self.end - self.start


def parse_path_info(path_name: str) -> tuple[str, str, str, str, str] | None:
    match = re.search(r"^(PAN\d+)#([12])#.*_(chr(?:[0-9]+|X|Y))_([pq]arm)$", path_name)
    if not match:
        return None
    sample, hap, chrom, arm = match.groups()
    chrarm = chrom + ("p" if arm == "parm" else "q")
    return sample, hap, chrom, arm, chrarm


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def load_communities(path: Path) -> dict[str, str]:
    communities: dict[str, str] = {}
    if not path.exists():
        return communities
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            arm = row.get("arm", "").replace("_", "")
            if arm.startswith("chr"):
                communities[arm] = row.get("community", "unknown") or "unknown"
    return communities


def fmt_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return "NA"
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


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


def proportion(n: int, d: int) -> str:
    return f"{n / d:.6f}" if d else "NA"


def in_range(value: int, bounds: tuple[int, int]) -> bool:
    return bounds[0] <= value <= bounds[1]


def classify_hits(hits: Iterable[Hit]) -> str:
    hits = list(hits)
    if not hits:
        return "unresolved_no_call"
    hap_keys = {hit.hap_key for hit in hits}
    arms = {hit.arm_key for hit in hits}
    communities = {hit.community for hit in hits if hit.community != "unknown"}
    if len(hap_keys) == 1:
        return "unique_donor_haplotype"
    if len(arms) == 1:
        return "unique_donor_arm_multi_haplotype"
    if len(communities) == 1:
        return "same_community_ambiguous"
    return "cross_community_ambiguous"


def class_rank(resolvability: str) -> int:
    return {
        "unique_donor_haplotype": 0,
        "unique_donor_arm_multi_haplotype": 1,
        "same_community_ambiguous": 2,
        "cross_community_ambiguous": 3,
        "unresolved_no_call": 4,
    }.get(resolvability, 4)


def tract_resolvability(segments: list[Segment]) -> str:
    return classify_hits(hit for segment in segments for hit in segment.hits)


def compatibility(left: Segment, right: Segment, mode: str) -> tuple[bool, str]:
    if mode == "none":
        return left.end == right.start, "touching_only"
    if mode == "exact":
        shared = left.hap_keys & right.hap_keys
        return bool(shared), "shared_donor_haplotype" if shared else "incompatible_donor_haplotype"
    if mode == "arm":
        shared_haps = left.hap_keys & right.hap_keys
        if shared_haps:
            return True, "shared_donor_haplotype"
        shared_arms = left.arm_keys & right.arm_keys
        return bool(shared_arms), "shared_donor_arm" if shared_arms else "incompatible_donor_arm"
    if mode == "community":
        shared_haps = left.hap_keys & right.hap_keys
        if shared_haps:
            return True, "shared_donor_haplotype"
        shared_arms = left.arm_keys & right.arm_keys
        if shared_arms:
            return True, "shared_donor_arm"
        shared_communities = left.communities & right.communities
        return (
            bool(shared_communities),
            "shared_donor_community" if shared_communities else "incompatible_donor_community",
        )
    raise ValueError(f"unknown bridge mode: {mode}")


def load_segments(
    bed_path: Path,
    dataset: str,
    setting: str,
    pair: str,
    label: str,
    query_prefix: str,
    communities: dict[str, str],
    top_n: int,
    score_delta: float,
    min_score: float,
) -> list[Segment]:
    grouped: dict[tuple[str, int, int], list[tuple[float, int, Hit]]] = defaultdict(list)
    with open_text(bed_path) as handle:
        header = next(handle, "")
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or not fields[0].startswith(query_prefix):
                continue
            query_info = parse_path_info(fields[0])
            ref_info = parse_path_info(fields[3])
            if query_info is None or ref_info is None:
                continue
            score = float(fields[6])
            nth_best = int(float(fields[9]))
            sample, hap, ref_chrom, _, ref_chrarm = ref_info
            hit = Hit(
                sample=sample,
                hap=hap,
                chrom=ref_chrom,
                chrarm=ref_chrarm,
                community=communities.get(ref_chrarm, "unknown"),
                score=score,
                nth_best=nth_best,
                ref_name=fields[3],
            )
            key = (fields[0], int(fields[1]), int(fields[2]))
            grouped[key].append((score, nth_best, hit))

    segments: list[Segment] = []
    for (query_name, start, end), rows in grouped.items():
        query_info = parse_path_info(query_name)
        if query_info is None:
            continue
        _, _, query_chr, _, query_chrarm = query_info
        best = max(score for score, _, _ in rows)
        keep = [
            hit
            for score, nth_best, hit in rows
            if nth_best <= top_n and score >= min_score and (best - score) <= score_delta
        ]
        if not keep:
            continue
        keep.sort(key=lambda hit: (hit.nth_best, -hit.score, hit.sample, hit.hap, hit.chrarm))
        segments.append(
            Segment(
                dataset=dataset,
                setting=setting,
                pair=pair,
                label=label,
                query_name=query_name,
                query_chr=query_chr,
                query_chrarm=query_chrarm,
                start=start,
                end=end,
                best_score=best,
                hits=keep,
            )
        )
    segments.sort(key=lambda segment: (segment.query_name, segment.start, segment.end))
    return segments


def call_first_best_tracts(segments: list[Segment]) -> list[Tract]:
    n1_segments: list[Segment] = []
    for segment in segments:
        first = min(segment.hits, key=lambda hit: (hit.nth_best, -hit.score, hit.sample, hit.hap, hit.chrarm))
        n1_segments.append(
            Segment(
                dataset=segment.dataset,
                setting=segment.setting,
                pair=segment.pair,
                label=segment.label,
                query_name=segment.query_name,
                query_chr=segment.query_chr,
                query_chrarm=segment.query_chrarm,
                start=segment.start,
                end=segment.end,
                best_score=segment.best_score,
                hits=[first],
            )
        )
    return merge_segments(n1_segments, max_bridge_gap=0, bridge_mode="exact")


def merge_segments(segments: list[Segment], max_bridge_gap: int, bridge_mode: str) -> list[Tract]:
    tracts: list[Tract] = []
    for query_name, query_segments in groupby_query(segments).items():
        current: Tract | None = None
        for segment in query_segments:
            if current is None:
                current = start_tract(segment)
                continue
            gap = segment.start - current.end
            compatible, reason = compatibility(current.segments[-1], segment, bridge_mode)
            if gap >= 0 and gap <= max_bridge_gap and compatible:
                current.end = max(current.end, segment.end)
                current.segments.append(segment)
                current.bridge_bp += gap
                current.bridge_reasons[reason] += 1
            else:
                tracts.append(current)
                current = start_tract(segment)
        if current is not None:
            tracts.append(current)
    return tracts


def start_tract(segment: Segment) -> Tract:
    return Tract(
        dataset=segment.dataset,
        setting=segment.setting,
        pair=segment.pair,
        label=segment.label,
        query_name=segment.query_name,
        query_chr=segment.query_chr,
        query_chrarm=segment.query_chrarm,
        start=segment.start,
        end=segment.end,
        segments=[segment],
        bridge_bp=0,
        bridge_reasons=Counter(),
    )


def groupby_query(segments: list[Segment]) -> dict[str, list[Segment]]:
    grouped: dict[str, list[Segment]] = defaultdict(list)
    for segment in segments:
        grouped[segment.query_name].append(segment)
    for values in grouped.values():
        values.sort(key=lambda segment: (segment.start, segment.end))
    return grouped


def tract_donor_metadata(tract: Tract) -> dict[str, str]:
    hits = [hit for segment in tract.segments for hit in segment.hits]
    hap_keys = sorted({hit.hap_key for hit in hits})
    arms = sorted({hit.arm_key for hit in hits})
    communities = sorted({hit.community for hit in hits})
    scores = [hit.score for hit in hits]
    return {
        "donor_haplotypes": ",".join(hap_keys),
        "donor_arms": ",".join(arms),
        "donor_communities": ",".join(communities),
        "min_score": fmt_float(min(scores) if scores else None),
        "max_score": fmt_float(max(scores) if scores else None),
    }


def interchromosomal(tract: Tract) -> bool:
    return any(hit.chrom != tract.query_chr for segment in tract.segments for hit in segment.hits)


def write_tracts(tracts: list[Tract], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset",
        "setting",
        "pair",
        "label",
        "query_name",
        "query_chrarm",
        "start",
        "end",
        "length_bp",
        "n_segments",
        "resolvability",
        "is_interchromosomal",
        "bridge_bp",
        "bridge_reasons",
        "donor_haplotypes",
        "donor_arms",
        "donor_communities",
        "min_score",
        "max_score",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for tract in tracts:
            meta = tract_donor_metadata(tract)
            row = {
                "dataset": tract.dataset,
                "setting": tract.setting,
                "pair": tract.pair,
                "label": tract.label,
                "query_name": tract.query_name,
                "query_chrarm": tract.query_chrarm,
                "start": tract.start,
                "end": tract.end,
                "length_bp": tract.length,
                "n_segments": len(tract.segments),
                "resolvability": tract_resolvability(tract.segments),
                "is_interchromosomal": interchromosomal(tract),
                "bridge_bp": tract.bridge_bp,
                "bridge_reasons": ",".join(f"{k}:{v}" for k, v in sorted(tract.bridge_reasons.items())),
            }
            row.update(meta)
            writer.writerow(row)


def summary_row(
    dataset: str,
    setting: str,
    tracts: list[Tract],
    baseline_first_best_n: int,
    high_conf_n: int | None,
    lower_n1_n: int | None,
) -> dict[str, str]:
    inter = [tract for tract in tracts if interchromosomal(tract)]
    lengths = [tract.length for tract in inter]
    n = len(lengths)
    resolvability = Counter(tract_resolvability(tract.segments) for tract in inter)
    bridged = sum(1 for tract in inter if tract.bridge_bp > 0)
    return {
        "dataset": dataset,
        "setting": setting,
        "n_tracts": str(n),
        "baseline_first_best_n": str(baseline_first_best_n),
        "high_conf_m1000_patch_n": str(high_conf_n) if high_conf_n is not None else "NA",
        "lower_m0_n1_firstbest_n": str(lower_n1_n) if lower_n1_n is not None else "NA",
        "merged_vs_first_best_n": str(max(0, baseline_first_best_n - n)),
        "split_vs_first_best_n": str(max(0, n - baseline_first_best_n)),
        "ambiguous_tract_n": str(
            resolvability["unique_donor_arm_multi_haplotype"]
            + resolvability["same_community_ambiguous"]
            + resolvability["cross_community_ambiguous"]
            + resolvability["unresolved_no_call"]
        ),
        "unique_donor_haplotype_n": str(resolvability["unique_donor_haplotype"]),
        "unique_donor_arm_multi_haplotype_n": str(resolvability["unique_donor_arm_multi_haplotype"]),
        "same_community_ambiguous_n": str(resolvability["same_community_ambiguous"]),
        "cross_community_ambiguous_n": str(resolvability["cross_community_ambiguous"]),
        "unresolved_no_call_n": str(resolvability["unresolved_no_call"]),
        "bridged_tract_n": str(bridged),
        "min_bp": fmt_float(float(min(lengths)) if lengths else None),
        "q25_bp": fmt_float(quantile(lengths, 0.25)),
        "median_bp": fmt_float(quantile(lengths, 0.50)),
        "q75_bp": fmt_float(quantile(lengths, 0.75)),
        "iqr_bp": fmt_float(
            (quantile(lengths, 0.75) or math.nan) - (quantile(lengths, 0.25) or math.nan)
            if lengths
            else None
        ),
        "max_bp": fmt_float(float(max(lengths)) if lengths else None),
        "conversion_22_95_n": str(sum(in_range(length, CONVERSION_RANGE) for length in lengths)),
        "conversion_22_95_prop": proportion(sum(in_range(length, CONVERSION_RANGE) for length in lengths), n),
        "near_conversion_11_190_n": str(sum(in_range(length, NEAR_CONVERSION_RANGE) for length in lengths)),
        "near_conversion_11_190_prop": proportion(
            sum(in_range(length, NEAR_CONVERSION_RANGE) for length in lengths), n
        ),
        "crossover_318_688_n": str(sum(in_range(length, CROSSOVER_RANGE) for length in lengths)),
        "crossover_318_688_prop": proportion(sum(in_range(length, CROSSOVER_RANGE) for length in lengths), n),
        "near_crossover_159_1376_n": str(sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths)),
        "near_crossover_159_1376_prop": proportion(
            sum(in_range(length, NEAR_CROSSOVER_RANGE) for length in lengths), n
        ),
    }


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_overall_n(path: Path, group_type: str) -> int | None:
    if not path.exists():
        return None
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("group_type") == group_type:
                return int(row["n"])
    return None


def parse_run(spec: str) -> dict[str, object]:
    parts = spec.split(":")
    if len(parts) != 7:
        raise ValueError(
            "run specs must be name:merge_m:top_n:score_delta:min_score:max_bridge_gap:bridge_mode"
        )
    name, merge_m, top_n, score_delta, min_score, max_bridge_gap, bridge_mode = parts
    return {
        "name": name,
        "merge_m": merge_m,
        "top_n": int(top_n),
        "score_delta": float(score_delta),
        "min_score": float(min_score),
        "max_bridge_gap": int(max_bridge_gap),
        "bridge_mode": bridge_mode,
    }


def bed_for_pair(untangle_dir: Path, pair: str, merge_m: str) -> Path:
    return untangle_dir / f"{pair}.e50000.m{merge_m}.bed.gz"


def write_visual_inputs(tracts: list[Tract], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "representative_segments.tsv"
    by_class: dict[str, list[Tract]] = defaultdict(list)
    for tract in tracts:
        if interchromosomal(tract):
            by_class[tract_resolvability(tract.segments)].append(tract)
    selected: list[Tract] = []
    for resolvability in [
        "unique_donor_haplotype",
        "unique_donor_arm_multi_haplotype",
        "same_community_ambiguous",
        "cross_community_ambiguous",
        "unresolved_no_call",
    ]:
        ranked = sorted(by_class.get(resolvability, []), key=lambda tract: (-tract.bridge_bp, -tract.length))
        selected.extend(ranked[:3])
    selected = selected[:15]
    fieldnames = [
        "setting",
        "pair",
        "query_name",
        "tract_start",
        "tract_end",
        "segment_start",
        "segment_end",
        "resolvability",
        "segment_resolvability",
        "bridge_bp",
        "donor_arms",
        "donor_communities",
        "best_score",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for tract in selected:
            meta = tract_donor_metadata(tract)
            for segment in tract.segments:
                writer.writerow(
                    {
                        "setting": tract.setting,
                        "pair": tract.pair,
                        "query_name": tract.query_name,
                        "tract_start": tract.start,
                        "tract_end": tract.end,
                        "segment_start": segment.start,
                        "segment_end": segment.end,
                        "resolvability": tract_resolvability(tract.segments),
                        "segment_resolvability": classify_hits(segment.hits),
                        "bridge_bp": tract.bridge_bp,
                        "donor_arms": meta["donor_arms"],
                        "donor_communities": meta["donor_communities"],
                        "best_score": fmt_float(segment.best_score),
                    }
                )
    return path


def write_plot_script(out_dir: Path) -> Path:
    script = out_dir / "plot_representative_segments.py"
    script.write_text(
        """#!/usr/bin/env python3
import csv
import html
from pathlib import Path

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "representative_segments.tsv").open(), delimiter="\\t"))
if not rows:
    raise SystemExit("representative_segments.tsv is empty")

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
left = 90
right = 500
row_h = 28
top = 52
height = top + row_h * len(keys) + 70
min_x = min(int(row["tract_start"]) for row in rows)
max_x = max(int(row["tract_end"]) for row in rows)
span = max(1, max_x - min_x)

def sx(value):
    return left + (int(value) - min_x) / span * (width - left - right)

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
    "<style>text{font-family:Arial,sans-serif;font-size:11px}.title{font-size:16px;font-weight:bold}.axis{fill:#333}.label{font-size:9px}</style>",
    f'<text class="title" x="{left}" y="24">Representative multimap-aware WashU untangle tract calls</text>',
    f'<line x1="{left}" x2="{width-right}" y1="{height-42}" y2="{height-42}" stroke="#333" stroke-width="1"/>',
    f'<text class="axis" x="{left}" y="{height-20}">child/query coordinate within 500 kb flank (bp)</text>',
]
for y, key in enumerate(keys):
    group = [row for row in rows if (row["setting"], row["pair"], row["query_name"], row["tract_start"], row["tract_end"]) == key]
    tract_start = int(key[3])
    tract_end = int(key[4])
    cy = top + y * row_h
    parts.append(f'<line x1="{sx(tract_start):.1f}" x2="{sx(tract_end):.1f}" y1="{cy}" y2="{cy}" stroke="#222" stroke-width="1" opacity="0.5"/>')
    for row in group:
        start = int(row["segment_start"])
        end = int(row["segment_end"])
        cls = row["segment_resolvability"]
        parts.append(
            f'<rect x="{sx(start):.1f}" y="{cy-6}" width="{max(1, sx(end)-sx(start)):.1f}" height="12" '
            f'fill="{colors.get(cls, "#999999")}"><title>{html.escape(cls)}; {start}-{end}</title></rect>'
        )
    label = f"{key[0]} {key[1]}:{key[2].split('#')[-1]} {tract_end - tract_start} bp"
    parts.append(f'<text class="label" x="{width-right+12}" y="{cy+4}">{html.escape(label)}</text>')

legend_y = height - 60
legend_x = left
for cls, color in colors.items():
    parts.append(f'<rect x="{legend_x}" y="{legend_y}" width="12" height="12" fill="{color}"/>')
    parts.append(f'<text class="label" x="{legend_x+17}" y="{legend_y+10}">{html.escape(cls)}</text>')
    legend_x += 225
    if legend_x > width - 260:
        legend_x = left
        legend_y += 18
parts.append("</svg>")
(HERE / "representative_segments.svg").write_text("\\n".join(parts))
"""
    )
    script.chmod(0o755)
    return script


def write_markdown(
    rows: list[dict[str, str]],
    tracts_path: Path,
    visual_tsv: Path,
    plot_script: Path,
    out_md: Path,
    sweepga_help: str,
) -> None:
    default = rows[0]
    bridge = rows[1] if len(rows) > 1 else rows[0]
    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(REPO_ROOT))
        except ValueError:
            return str(path)

    visual_svg = visual_tsv.with_suffix(".svg")
    text = f"""# WashU Untangle Multimap-Aware Tract Calls

Script: `scripts/pedigree/untangle_multimap_tracts.py`

Primary outputs:

- `scripts/pedigree/untangle_multimap_tracts.tsv`
- `scripts/pedigree/untangle_multimap_tract_summary.tsv`
- `{rel(visual_tsv)}`
- `{rel(visual_svg)}`
- `{rel(plot_script)}`

## What Was Tried

The previous lower-threshold check merged only `nth.best == 1` intervals by a
single donor haplotype. This pass instead groups every exact child/query
interval in the WashU `odgi untangle` BEDs, keeps top-N hits that are within a
score delta of the best hit, and treats those hits as the interval's donor
equivalence class. Adjacent intervals are merged only when those equivalence
classes remain compatible under the selected bridge mode.

The default sensitivity runs are:

| setting | top-N | score delta | min score | max bridge gap | bridge mode |
|---|---:|---:|---:|---:|---|
"""
    for row in rows:
        setting = row["setting"]
        if "bridge1kb" in setting:
            text += f"| {setting} | 8 | 0.002 | 0.8 | 1000 | community |\n"
        else:
            text += f"| {setting} | 4 | 0.001 | 0.8 | 0 | arm |\n"
    text += f"""
The committed run processes the available m1000 BEDs directly. The script also
accepts lower-threshold BEDs through the `merge_m` field in `--run` specs, for
example `--run lower_m0:0:8:0.002:0.8:1000:community`, but those large m0 BED
intermediates are not committed in this worktree. For that reason the summary
table carries forward the existing `m0/n1` first-best lower-merge denominator
from `scripts/pedigree/patch_tract_lower_merge_summary.tsv` as a comparison
rather than silently treating it as multimap-aware evidence.

## sweepga Check

`sweepga --help` shows that sweepga accepts FASTA inputs or a single PAF and
then applies plane-sweep filtering with `--num-mappings`. The WashU evidence
available here is already the downstream `odgi untangle` BED output, with
child/query coordinates, reference path coordinates, scores, and `nth.best`.
There is no PAF-equivalent retained in the repository or the WashU untangle
directory. Re-running sweepga would therefore require reconstructing a separate
FASTA/PAF alignment path rather than preserving odgi's graph-specific interval
calls. For this task, the appropriate equivalent is an interval sweep over the
untangle BED rows themselves.

Relevant sweepga interface excerpt:

```text
{sweepga_help.strip()}
```

## Results

Under the recommended default `{default['setting']}`, the multimap-aware caller
returns {default['n_tracts']} interchromosomal candidate tracts. The matched
first-best projection contains {default['baseline_first_best_n']} tracts, so
the multimap-aware representation merges {default['merged_vs_first_best_n']} and
splits {default['split_vs_first_best_n']} relative to the arbitrary first-best
view. Ambiguity is explicit: {default['ambiguous_tract_n']} tracts are not
unique donor haplotypes, including {default['unique_donor_arm_multi_haplotype_n']}
unique-arm/multiple-haplotype tracts, {default['same_community_ambiguous_n']}
same-community ambiguous tracts, and
{default['cross_community_ambiguous_n']} cross-community ambiguous tracts.

With the more permissive `{bridge['setting']}` sensitivity, the caller returns
{bridge['n_tracts']} interchromosomal candidate tracts, merges
{bridge['merged_vs_first_best_n']} relative to first-best, and records
{bridge['bridged_tract_n']} tracts with non-zero bridged sequence. This confirms
that consecutive m1000 intervals can be merged through compatible donor
equivalence classes. In this m1000 run the compatible joins are mostly zero-gap
adjacency rather than non-zero missing sequence; some intervals remain genuinely
unresolved or cross-community ambiguous.

Length distributions against the primate literature windows are summarized in
`scripts/pedigree/untangle_multimap_tract_summary.tsv`. For the recommended
default, median tract length is {default['median_bp']} bp (IQR
{default['iqr_bp']} bp; min {default['min_bp']} bp; max {default['max_bp']} bp).
Counts in the descriptive windows are:

- 22-95 bp: {default['conversion_22_95_n']}/{default['n_tracts']}
  ({default['conversion_22_95_prop']})
- 318-688 bp: {default['crossover_318_688_n']}/{default['n_tracts']}
  ({default['crossover_318_688_prop']})
- 159-1376 bp: {default['near_crossover_159_1376_n']}/{default['n_tracts']}
  ({default['near_crossover_159_1376_prop']})

## Interpretation

The analysis no longer treats `nth.best == 1` as uniquely true when tied or
near-tied donor paths exist. It strengthens the tract-length audit by showing
which m1000 intervals are mergeable under explicit donor equivalence classes,
and by identifying same-arm, same-community, and cross-community ambiguity.

The result remains a supportive compatibility analysis. It does not by itself
prove conversion or crossover mechanisms, because repeated subtelomeric
haplotypes often make the exact donor unresolved. The useful claim is cautious:
some WashU candidate regions are compatible with merged tract interpretations
once multimapping is represented, while a substantial fraction should remain
ambiguous rather than be forced into a first-best donor.
"""
    out_md.write_text(text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--untangle-dir", type=Path, default=UNTANGLE_DIR)
    parser.add_argument("--communities", type=Path, default=COMMUNITY_TSV)
    parser.add_argument("--out-tracts", type=Path, default=OUT_TRACTS)
    parser.add_argument("--out-summary", type=Path, default=OUT_SUMMARY)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--out-md", type=Path, default=OUT_MD)
    parser.add_argument(
        "--run",
        action="append",
        default=list(DEFAULT_RUNS),
        help="name:merge_m:top_n:score_delta:min_score:max_bridge_gap:bridge_mode",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    communities = load_communities(args.communities)
    high_conf_n = read_overall_n(PATCH_SUMMARY_TSV, "overall")
    lower_n1_n = read_overall_n(LOWER_SUMMARY_TSV, "lower_merge_overall")
    all_tracts: list[Tract] = []
    summary_rows: list[dict[str, str]] = []
    for run_spec in args.run:
        run = parse_run(run_spec)
        setting = str(run["name"])
        run_segments: list[Segment] = []
        for pair, (query_prefix, label) in PAIRINGS.items():
            bed_path = bed_for_pair(args.untangle_dir, pair, str(run["merge_m"]))
            if not bed_path.exists():
                print(f"Skipping missing BED: {bed_path}")
                continue
            run_segments.extend(
                load_segments(
                    bed_path=bed_path,
                    dataset=f"e50000.m{run['merge_m']}",
                    setting=setting,
                    pair=pair,
                    label=label,
                    query_prefix=query_prefix,
                    communities=communities,
                    top_n=int(run["top_n"]),
                    score_delta=float(run["score_delta"]),
                    min_score=float(run["min_score"]),
                )
            )
        first_best = call_first_best_tracts(run_segments)
        tracts = merge_segments(
            run_segments,
            max_bridge_gap=int(run["max_bridge_gap"]),
            bridge_mode=str(run["bridge_mode"]),
        )
        all_tracts.extend(tracts)
        summary_rows.append(
            summary_row(
                dataset=f"e50000.m{run['merge_m']}",
                setting=setting,
                tracts=tracts,
                baseline_first_best_n=sum(1 for tract in first_best if interchromosomal(tract)),
                high_conf_n=high_conf_n,
                lower_n1_n=lower_n1_n,
            )
        )
    if not all_tracts:
        raise RuntimeError("No tracts were called; check input BED paths and run specs")
    write_tracts(all_tracts, args.out_tracts)
    write_tsv(summary_rows, args.out_summary)
    visual_tsv = write_visual_inputs(all_tracts, args.out_dir)
    plot_script = write_plot_script(args.out_dir)
    try:
        subprocess.run(["python3", str(plot_script)], check=True)
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        print(f"Plot generation skipped/failed: {exc}")
    sweepga = subprocess.run(
        ["/home/erikg/.cargo/bin/sweepga", "--help"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    sweepga_excerpt = "\n".join(
        line
        for line in sweepga.stdout.splitlines()
        if "--num-mappings" in line or "Input files" in line or "FASTA" in line or "PAF" in line
    )
    write_markdown(summary_rows, args.out_tracts, visual_tsv, plot_script, args.out_md, sweepga_excerpt)
    print(f"Wrote {args.out_tracts}")
    print(f"Wrote {args.out_summary}")
    print(f"Wrote {args.out_md}")
    print(f"Wrote {visual_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
