#!/usr/bin/env python3
import glob
import os
import sys

from common import PACKAGE_DIR, overlap_bp, paf_records, read_tsv, write_tsv


def summarize_one(path):
    records = list(paf_records(path))
    total_bp = sum(r["query_end"] - r["query_start"] for r in records)
    inter = [r for r in records if r["query_arm"] != "unknown" and r["target_arm"] != "unknown" and r["query_arm"] != r["target_arm"]]
    best_by_query_arm = {}
    for r in records:
        key = (r["query_name"], r["target_hap"], r["target_arm"])
        bp = r["query_end"] - r["query_start"]
        best_by_query_arm[key] = best_by_query_arm.get(key, 0) + bp
    top = sorted(best_by_query_arm.items(), key=lambda kv: kv[1], reverse=True)[:12]
    return {
        "file": path,
        "records": str(len(records)),
        "query_aligned_bp_sum": str(total_bp),
        "inter_arm_records": str(len(inter)),
        "inter_arm_query_bp_sum": str(sum(r["query_end"] - r["query_start"] for r in inter)),
        "top_query_target_arm_bp": ";".join("%s|%s|%s:%d" % (k[0], k[1], k[2], v) for k, v in top),
    }


def write_segment_support(paf_paths):
    selected = read_tsv(os.path.join(PACKAGE_DIR, "..", "fig5_synteny_recombination_schematic", "selected_segments.tsv"))
    comparison_to_files = {}
    for path in paf_paths:
        base = os.path.basename(path)
        cid = base.split(".", 1)[0]
        comparison_to_files.setdefault(cid, []).append(path)
    comp_rows = read_tsv(os.path.join(PACKAGE_DIR, "config", "comparisons.tsv"))
    pair_to_comparisons = {}
    for row in comp_rows:
        pair_to_comparisons.setdefault(row["pair"], []).append(row["comparison_id"])

    rows = []
    for seg in selected:
        expected_q = seg["query_name"]
        expected_t = seg["target_name"]
        q_start = int(seg["local_query_start_0based"])
        q_end = int(seg.get("local_query_end_0based_exclusive") or seg.get("local_query_end_0based"))
        expected_target_arm = seg["target_arm"]
        expected_target_hap = seg["target_haplotype"]
        best = None
        candidate_files = []
        for cid in pair_to_comparisons.get(seg["pair"], []):
            candidate_files += comparison_to_files.get(cid, [])
        for path in candidate_files:
            filter_id = os.path.basename(path).replace(path.split(os.sep)[-1].split(".")[0] + ".", "").replace(".paf.gz", "")
            for r in paf_records(path):
                if r["query_name"] != expected_q:
                    continue
                ov = overlap_bp(q_start, q_end, r["query_start"], r["query_end"])
                if ov == 0:
                    continue
                same_hap = r["target_hap"] == expected_target_hap
                same_arm = r["target_arm"] == expected_target_arm
                score = (1 if same_hap else 0, 1 if same_arm else 0, ov, r["identity"], r["aln_length"])
                if best is None or score > best[0]:
                    best = (score, path, filter_id, r, ov)
        if best:
            _, path, filter_id, r, ov = best
            support = "agree" if r["target_hap"] == expected_target_hap and r["target_arm"] == expected_target_arm else "discordant_target"
            if ov < max(1, int(0.5 * (q_end - q_start))):
                support = "inconclusive_partial_overlap" if support == "agree" else support
            if support == "agree" and r["identity"] >= 0.95:
                recommendation = "direct_sweepga_can_be_primary_for_segment"
            elif support == "discordant_target":
                recommendation = "keep_graph_primary_pending_manual_review"
            else:
                recommendation = "inconclusive_keep_graph_primary"
            rows.append({
                "event_id": seg["event_id"],
                "pair": seg["pair"],
                "event_role": seg["event_role"],
                "query_name": expected_q,
                "query_local_interval": "%s-%s" % (q_start, q_end),
                "expected_target_name": expected_t,
                "expected_target_haplotype": expected_target_hap,
                "expected_target_arm": expected_target_arm,
                "direct_support": support,
                "direct_file": path,
                "direct_filter_id": filter_id,
                "direct_target_name": r["target_name"],
                "direct_target_haplotype": r["target_hap"],
                "direct_target_arm": r["target_arm"],
                "direct_query_interval": "%d-%d" % (r["query_start"], r["query_end"]),
                "overlap_bp": str(ov),
                "identity": "%.6f" % r["identity"],
                "aln_length": str(r["aln_length"]),
                "evidence_source_recommendation": recommendation,
            })
        else:
            rows.append({
                "event_id": seg["event_id"],
                "pair": seg["pair"],
                "event_role": seg["event_role"],
                "query_name": expected_q,
                "query_local_interval": "%s-%s" % (q_start, q_end),
                "expected_target_name": expected_t,
                "expected_target_haplotype": expected_target_hap,
                "expected_target_arm": expected_target_arm,
                "direct_support": "no_direct_overlap",
                "direct_file": "",
                "direct_filter_id": "",
                "direct_target_name": "",
                "direct_target_haplotype": "",
                "direct_target_arm": "",
                "direct_query_interval": "",
                "overlap_bp": "0",
                "identity": "",
                "aln_length": "",
                "evidence_source_recommendation": "no_direct_signal_keep_graph_primary",
            })
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "direct_vs_graph_concordance.tsv"), rows,
              ["event_id", "pair", "event_role", "query_name", "query_local_interval",
               "expected_target_name", "expected_target_haplotype", "expected_target_arm",
               "direct_support", "direct_file", "direct_filter_id", "direct_target_name",
               "direct_target_haplotype", "direct_target_arm", "direct_query_interval",
               "overlap_bp", "identity", "aln_length", "evidence_source_recommendation"])


def main():
    paf_paths = sorted(glob.glob(os.path.join(PACKAGE_DIR, "raw_paf", "*.paf.gz")) +
                       glob.glob(os.path.join(PACKAGE_DIR, "filtered_paf", "*.paf.gz")))
    rows = [summarize_one(path) for path in paf_paths]
    write_tsv(os.path.join(PACKAGE_DIR, "summaries", "paf_file_summary.tsv"), rows,
              ["file", "records", "query_aligned_bp_sum", "inter_arm_records", "inter_arm_query_bp_sum", "top_query_target_arm_bp"])
    write_segment_support(paf_paths)
    return 0


if __name__ == "__main__":
    sys.exit(main())
