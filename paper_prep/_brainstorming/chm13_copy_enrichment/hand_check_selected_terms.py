#!/usr/bin/env python3
"""Direct coordinate hand checks for prespecified duplicated families/terms.

No engine code or cached target flags are used to select loci.  The script
joins raw half-open target intervals to physical-locus coordinates, then joins
the frozen engine collections by locus_id and writes the contributing IDs.
"""

import csv
import gzip
from pathlib import Path


HERE = Path(__file__).resolve().parent
CHECKS = (
    ("HGNC_group", "HGNC_GROUP:151"),       # OR4F family
    ("HGNC_group", "HGNC_GROUP:3661"),      # tubulin beta family
    ("HGNC_group", "HGNC_GROUP:521"),       # PRD homeoboxes/pseudogenes
    ("GO_MF", "GO:0004984"),                # olfactory receptor activity
    ("biotype", "BIOTYPE:pseudogene"),
    ("biotype", "BIOTYPE:protein_coding"),
)


def read(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def selected(loci, intervals, assignment):
    result = set()
    for locus_id, locus in loci.items():
        for interval in intervals:
            if locus["chromosome"] != interval["chromosome"] or locus["arm"] != interval["arm"]:
                continue
            if assignment == "midpoint":
                hit = int(interval["start0"]) <= int(locus["midpoint0"]) < int(interval["end0"])
            else:
                hit = int(locus["start0"]) < int(interval["end0"]) and int(interval["start0"]) < int(locus["end0"])
            if hit:
                result.add(locus_id)
                break
    return result


def main():
    loci = {row["locus_id"]: row for row in read(HERE / "analysis_ready/chm13_gene_loci.tsv.gz")}
    intervals = list(read(HERE / "analysis_ready/chm13_phr_intervals.tsv"))
    midpoint = selected(loci, intervals, "midpoint")
    overlap = selected(loci, intervals, "overlap")
    if len(intervals) != 37 or len(midpoint) != 402 or len(overlap) != 412:
        raise AssertionError("unexpected target template/locus counts")
    if midpoint != {key for key, row in loci.items() if row["in_phr_midpoint"] == "1"}:
        raise AssertionError("direct midpoint recount disagrees with prepared target flags")
    if overlap != {key for key, row in loci.items() if row["in_phr_any_overlap"] == "1"}:
        raise AssertionError("direct overlap recount disagrees with prepared target flags")

    rows_out = []
    for collection, term_id in CHECKS:
        term_rows = [row for row in read(HERE / "engine_terms" / (collection + ".tsv.gz"))
                     if row["term_id"] == term_id]
        if not term_rows:
            raise AssertionError("prespecified term missing: %s/%s" % (collection, term_id))
        all_ids = {row["locus_id"] for row in term_rows}
        for assignment, chosen in (("midpoint", midpoint), ("overlap", overlap)):
            ids = sorted(all_ids & chosen)
            rows_out.append({
                "collection": collection, "term_id": term_id,
                "term_name": term_rows[0]["term_name"], "assignment": assignment,
                "genome_copy_count": len(all_ids), "observed_copy_count": len(ids),
                "observed_arm_breadth": len({loci[item]["arm"] for item in ids}),
                "physical_locus_ids": ",".join(ids),
                "copy_ids": ",".join(sorted(row["copy_id"] for row in term_rows
                                              if row["locus_id"] in chosen)),
                "frozen_source": term_rows[0]["frozen_source"],
            })
    output = HERE / "hand_check_selected_terms.tsv"
    fields = list(rows_out[0])
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)
    print("wrote %d checks; 37 intervals, %d midpoint loci, %d overlap loci" %
          (len(rows_out), len(midpoint), len(overlap)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
