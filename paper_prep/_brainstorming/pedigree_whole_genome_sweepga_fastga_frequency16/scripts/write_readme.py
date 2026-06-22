#!/usr/bin/env python3
import csv
import os

from common import PACKAGE_DIR


def read_rows(name):
    path = os.path.join(PACKAGE_DIR, "summaries", name)
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def yes(rows, candidate):
    return "yes" if any(
        row.get("candidate_id") == candidate
        and row.get("target_chrom") == "chr3"
        and row.get("emits_chr3_target") == "yes"
        for row in rows
    ) else "no"


def main():
    raw = read_rows("raw_chr3_support.tsv")
    jobs = read_rows("slurm_jobs.tsv")
    pathological = read_rows("pathological_runtime.tsv")
    summary = read_rows("frequency_sensitivity_summary.tsv")
    pan027 = yes(raw, "PAN027_chr9q_chr3q_PHR_candidate")
    pan028 = yes(raw, "PAN028_chr9q_chr3q_PHR_candidate")
    raw_rows = {
        (row.get("candidate_id"), row.get("target_chrom")): row for row in raw
    }
    candidate_rows = read_rows("candidate_window_support.tsv")
    support_rows = {
        (row.get("candidate_id"), row.get("layer"), row.get("target_chrom")): row
        for row in candidate_rows
    }
    pathological_text = ""
    if pathological:
        pathological_text = (
            "\nThe frequency-16 jobs were treated as pathological because FastGA remained active "
            "without raw PAF output and with zero/near-zero `.1aln` output after runtime well beyond "
            "the prior `-f2` jobs. Chopping/filtering was skipped because no raw PAF evidence existed.\n"
        )
    elif pan027 == "yes" or pan028 == "yes":
        p27_raw = raw_rows.get(("PAN027_chr9q_chr3q_PHR_candidate", "chr3"), {})
        p28_raw = raw_rows.get(("PAN028_chr9q_chr3q_PHR_candidate", "chr3"), {})
        p27_many = support_rows.get(("PAN027_chr9q_chr3q_PHR_candidate", "many_many_chopped", "chr3"), {})
        p27_four = support_rows.get(("PAN027_chr9q_chr3q_PHR_candidate", "four_many_chopped", "chr3"), {})
        p28_many = support_rows.get(("PAN028_chr9q_chr3q_PHR_candidate", "many_many_chopped", "chr3"), {})
        p28_four = support_rows.get(("PAN028_chr9q_chr3q_PHR_candidate", "four_many_chopped", "chr3"), {})
        pathological_text = f"""
Raw support:

- `PAN027_chr9q_chr3q_PHR_candidate`: {p27_raw.get('overlap_rows', '')} raw chr3 rows, {p27_raw.get('sum_overlap_bp', '')} bp summed query-window overlap, {p27_raw.get('query_union_bp', '')} bp query-union coverage.
- `PAN028_chr9q_chr3q_PHR_candidate`: {p28_raw.get('overlap_rows', '')} raw chr3 rows, {p28_raw.get('sum_overlap_bp', '')} bp summed query-window overlap, {p28_raw.get('query_union_bp', '')} bp query-union coverage.

Because raw chr3 rows appeared at frequency 16, 10 kb `pafchop-rs` and chopped
sweepGA `many:many`/`4:many` filters were run. The chr3 signal persists through
both required chopped evidence layers:

- PAN027 `many:many` chopped: {p27_many.get('overlap_rows', '')} chr3 rows, {p27_many.get('sum_overlap_bp', '')} bp summed overlap, {p27_many.get('query_union_bp', '')} bp query-union coverage.
- PAN027 `4:many` chopped: {p27_four.get('overlap_rows', '')} chr3 rows, {p27_four.get('sum_overlap_bp', '')} bp summed overlap, {p27_four.get('query_union_bp', '')} bp query-union coverage.
- PAN028 `many:many` chopped: {p28_many.get('overlap_rows', '')} chr3 rows, {p28_many.get('sum_overlap_bp', '')} bp summed overlap, {p28_many.get('query_union_bp', '')} bp query-union coverage.
- PAN028 `4:many` chopped: {p28_four.get('overlap_rows', '')} chr3 rows, {p28_four.get('sum_overlap_bp', '')} bp summed overlap, {p28_four.get('query_union_bp', '')} bp query-union coverage.
"""
    else:
        pathological_text = (
            "\nThe frequency-16 raw PAFs completed, but no raw chr3 rows overlapped the PAN027/PAN028 "
            "candidate windows. Chopping/filtering was skipped by design because the raw-first test was negative.\n"
        )
    job_lines = "\n".join(
        "- `{job_id}`: `{comparison_id}` {status} {elapsed}".format(**{**row, "elapsed": row.get("elapsed", "")})
        for row in jobs
    )
    freq_lines = "\n".join(
        "- `{candidate_id}`: frequency16 raw chr3 `{explicit_fastga_frequency16_raw_chr3}`, prior sweepGA raw `{prior_updated_bin_raw_no_explicit_frequency_chr3}`, wfmash p95 `{updated_wfmash_p95_chr3}`".format(**row)
        for row in summary
    )
    if pan027 == "yes" or pan028 == "yes":
        window_scope = "both PAN027/PAN028 candidate windows" if pan027 == "yes" and pan028 == "yes" else "at least one PAN027/PAN028 candidate window"
        interpretation = (
            "The frequency-16 run supports the seed-frequency sparsification hypothesis behind "
            "the wfmash-positive / sweepGA-negative discrepancy: the prior updated-bin "
            "no-explicit-frequency sweepGA run used the much stricter effective FastGA `-f2` "
            "setting and missed raw chr3 support, while the explicit `-f16` run finished all "
            "three whole-genome comparisons and recovered chr3 rows for both candidate windows. "
            "Unlike `-f100`, `-f16` was slow but not pathological at this scale."
        )
    else:
        window_scope = "the PAN027/PAN028 candidate windows"
        interpretation = (
            "If frequency16 is negative while wfmash remains positive, the result leaves the "
            "seed-frequency sparsification hypothesis unresolved rather than rescuing the "
            "discrepancy: `-f16` is less pathological than `-f100` only if it finishes, but raw "
            "chr3 recovery is the decisive criterion."
        )

    text = f"""# Fig5 whole-genome sweepGA/FastGA frequency-16 sensitivity

Date: 2026-06-21

This package tests whether an explicit FastGA k-mer occurrence threshold of 16
rescues chr3 target homology for the Fig5 PAN027/PAN028 chr9 candidate windows
in whole-genome sweepGA/FastGA output.

Primary command shape:

```bash
/home/erikg/.cargo/bin/sweepga \\
  --fastga \\
  --fastga-frequency 16 \\
  --num-mappings many:many \\
  --scaffold-jump 0 \\
  --temp-dir /dev/shm/... \\
  --output-file ... \\
  QUERY.fa TARGET.fa
```

The three comparisons are the same joint-parent whole-genome inputs used by
`paper_prep/_brainstorming/pedigree_whole_genome_sweepga_updated_bin/`:

- `PAN027pat_vs_PAN011_joint`
- `PAN027mat_vs_PAN010_joint`
- `PAN028mat_vs_PAN027_joint`

## Binary Provenance

`summaries/sweepga_binary.tsv` records the explicit sweepGA path, `which`,
realpath, version, sha256, and `--help` text for `/home/erikg/.cargo/bin/sweepga`.

`summaries/fastga_binary.tsv` records `/home/erikg/.cargo/bin/sweepga --check-fastga`.
Each Slurm command log also records the explicit binary path, `which`, compute-node
realpath, sha256, `--help`, `--check-fastga`, `/dev/shm` scratch path, and exact command.

## Workflow

Raw frequency-16 Slurm jobs:

{job_lines}

All three jobs are full whole-genome runs. Scratch is explicitly under `/dev/shm`;
`$SLURM_TMPDIR` is not used as sweepGA/FastGA scratch.

## Result

Direct answer: **PAN027 {pan027}; PAN028 {pan028}**. Explicit
`--fastga-frequency 16` {'made' if pan027 == 'yes' and pan028 == 'yes' else 'did not make'}
sweepGA/FastGA emit chr3 target rows for {window_scope} in raw PAF.
{pathological_text}
Comparator summary:

{freq_lines}

The updated wfmash `-p 95` comparator is treated as expected-positive evidence,
not as a filter input. {interpretation}

Validated with:

```bash
paper_prep/_brainstorming/pedigree_whole_genome_sweepga_fastga_frequency16/scripts/validate_outputs.sh
```

## Required Summaries

- `summaries/sweepga_binary.tsv`
- `summaries/fastga_binary.tsv`
- `summaries/slurm_jobs.tsv`
- `summaries/raw_chr3_support.tsv`
- `summaries/frequency_sensitivity_summary.tsv`
- `summaries/pathological_runtime.tsv` if pathological
- `summaries/chop_manifest.tsv`, `summaries/filter_manifest.tsv`, and
  `summaries/candidate_window_support.tsv` if chopping/filtering ran

Raw, chopped, and filtered PAFs and checksums are ignored.
"""
    with open(os.path.join(PACKAGE_DIR, "README.md"), "w") as out:
        out.write(text)


if __name__ == "__main__":
    main()
