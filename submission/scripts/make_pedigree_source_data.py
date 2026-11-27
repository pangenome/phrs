#!/usr/bin/env python3
"""Source Data for the pedigree Extended Data figures (one .xlsx per figure, one sheet per panel), from files in the repo.

usage (repo root): python3 submission/scripts/make_pedigree_source_data.py
ED4: every alignment record drawn in the dot plots (data/ed4_pairwise.paf; identity = matches / block length).
ED7 B: the 2 kb windows of panel B (make_fig5_zoom.R, SELF.p chr9 136.13-136.19 Mb): same-chromosome best match against the
  released assemblies (data/fig5_PAN027pat_vs_PAN011_joint...) and same-chromosome and non-homologous best matches against the
  corrected assemblies (data/fig5_corrected_v11_PAN027pat_vs_PAN011_joint...).
ED7 C: the 51 drawn non-homologous runs (scripts/pedigree/fig5_v131_check/results/fig5_final/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv);
  the grey same-chromosome chains come from the corrected class-winner table, which is in data/.
out: submission/source_data/Source_Data_ED4.xlsx, Source_Data_ED7.xlsx
"""
import csv, gzip, os
import openpyxl
OUT = "submission/source_data"; os.makedirs(OUT, exist_ok=True)
def save(name, sheets):
    wb = openpyxl.Workbook(); wb.remove(wb.active)
    for k, rows in sheets.items():
        ws = wb.create_sheet(k)
        for r in rows: ws.append(list(r))
    wb.save(f"{OUT}/{name}.xlsx")
    print(f"wrote {OUT}/{name}.xlsx: " + ", ".join(f"{k} {len(v) - 1} rows" for k, v in sheets.items()))
paf = [("query", "query_length", "query_start", "query_end", "strand", "target", "target_length", "target_start", "target_end", "matches", "block_length", "identity")]
for l in open("data/ed4_pairwise.paf"):
    f = l.split("\t"); paf.append((f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]), round(int(f[9]) / int(f[10]), 6)))
save("Source_Data_ED4", {"ED4": paf})
def win(path, release, classes):
    rows = []
    for r in csv.DictReader(gzip.open(path, "rt"), delimiter="\t"):
        if r["query_chrom"] == "chr9" and int(r["start"]) >= 136130000 and int(r["end"]) <= 136190000 and r["winner_class"] in classes:
            rows.append((release, int(r["start"]), int(r["end"]), r["winner_class"], r["other_seq"], float(r["estimated.identity"])))
    return rows
b = [("assemblies", "SELF.p_chr9_start", "SELF.p_chr9_end", "winner_class", "father_sequence", "best_match_identity")]
b += win("data/fig5_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz", "released", {"same_chrom"})
b += win("data/fig5_corrected_v11_PAN027pat_vs_PAN011_joint.class_winners.impg_similarity.tsv.gz", "corrected", {"same_chrom", "interchrom"})
c = [tuple(l.rstrip("\n").split("\t")) for l in open("scripts/pedigree/fig5_v131_check/results/fig5_final/PAN027pat_vs_PAN011_joint.whole_genome_ribbon_runs.tsv")]
save("Source_Data_ED7", {"ED7B": b, "ED7C": c})
