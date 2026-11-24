#!/usr/bin/env python3
"""Gene families with copies on several chromosome ends across the pangenome (R1-M1 test 1, family choice).

extract:   for each per-haplotype Liftoff GFF3 (PanSN contig names), the gene records that overlap a
           telomere-anchored 500 kb window (names contig:a-b_<arm>, a 1-based, from the windows FASTA .fai).
           Output one row per (window, gene): sample, hap, arm, gene_name, biotype, offset of the gene in the
           window, distance to the telomere end of the window, length, Liftoff coverage and identity.
summarize: per FAMILY, the arms that carry it and in how many haplotypes. A family is a gene_name, or, with
           --stem, the gene_name stripped of one copy suffix (see stem()). An arm carries
           a family when >= --min-frac of the haplotypes with a window on that arm have a copy in it.
"""
import argparse, csv, glob, gzip, os, re, statistics as st, sys, time
from collections import defaultdict

NAME = re.compile(r"^(.+):(\d+)-(\d+)_(chr[0-9XYM]+_(?:parm|qarm))$")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def stem(n):
    """the first rule that applies: olfactory subfamily (OR4F29 -> OR4F), P + number (SEPTIN14P22 -> SEPTIN14,
    IL9RP1 -> IL9R), number + optional letter + P after a letter (WASH7P -> WASH), number after L (DDX11L16 ->
    DDX11L), letter A-F after a digit (FAM138A -> FAM138, TUBB8B -> TUBB8); other names are kept whole"""
    for rx, rep in ((r"^(OR\d+[A-Z]+)\d+[A-Z]*P?$", r"\1"), (r"P\d+$", ""), (r"(?<=[A-Z])\d+[A-Z]?P$", ""),
                    (r"(?<=L)\d+$", ""), (r"(?<=\d)[A-F]$", "")):
        if re.search(rx, n):
            return re.sub(rx, rep, n)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("--fai", required=True, help="windows FASTA .fai")
    e.add_argument("--gff", nargs="+", required=True)
    e.add_argument("--out", required=True)
    s = sub.add_parser("summarize")
    s.add_argument("--tables", nargs="+", required=True)
    s.add_argument("--stem", action="store_true")
    s.add_argument("--min-frac", type=float, default=0.10)
    s.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "extract":
        win = defaultdict(list)
        for line in open(a.fai):
            n = line.split("\t")[0]
            m = NAME.match(n)
            if not m:
                raise SystemExit(f"FATAL: cannot parse window name {n!r}")
            win[m.group(1)].append((int(m.group(2)) - 1, int(m.group(3)), n, m.group(4)))
        log(f"{sum(map(len, win.values()))} windows on {len(win)} contigs")
        with open(a.out, "w") as fh:
            fh.write("window\tsample\thap\tarm\tgene_name\tbiotype\tstart\tend\ttelomere_distance\tlength\tcoverage\tidentity\n")
            for g in a.gff:
                t0, n_gene, n_kept, contigs = time.time(), 0, 0, set()
                with gzip.open(g, "rt") as gh:
                    for line in gh:
                        if "\tgene\t" not in line:
                            continue
                        f = line.rstrip("\n").split("\t")
                        if f[2] != "gene":
                            continue
                        n_gene += 1
                        contigs.add(f[0])
                        ws = win.get(f[0])
                        if not ws:
                            continue
                        gs, ge = int(f[3]) - 1, int(f[4])
                        for w0, w1, wn, arm in ws:
                            if gs < w1 and ge > w0:
                                at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
                                s0, s1 = max(gs, w0) - w0, min(ge, w1) - w0
                                L = w1 - w0
                                td = s0 if arm.endswith("parm") else L - s1
                                sm, hp = wn.split("#")[:2]
                                fh.write(f"{wn}\t{sm}\t{hp}\t{arm}\t{at.get('gene_name', at.get('ID'))}\t"
                                         f"{at.get('gene_biotype', 'NA')}\t{s0}\t{s1}\t{td}\t{ge - gs}\t"
                                         f"{at.get('coverage', 'NA')}\t{at.get('sequence_ID', 'NA')}\n")
                                n_kept += 1
                if not any(c in win for c in contigs):
                    raise SystemExit(f"FATAL: no contig of {g} matches a window contig (naming mismatch?)")
                log(f"{os.path.basename(g)}: {n_gene} genes, {n_kept} gene-window overlaps, {time.time() - t0:.1f} s")
        return
    rows = [r for t in a.tables for r in csv.DictReader(open(t), delimiter="\t")]
    log(f"{len(rows)} gene-window rows from {len(a.tables)} tables")
    hap_arm = defaultdict(set)          # arm -> haplotypes with a window that has >= 1 gene row
    fam = defaultdict(lambda: defaultdict(set))
    info = defaultdict(lambda: dict(names=set(), bt=defaultdict(int), lens=[], td=[]))
    for r in rows:
        k = stem(r["gene_name"]) if a.stem else r["gene_name"]
        h = (r["sample"], r["hap"])
        hap_arm[r["arm"]].add(h)
        fam[k][r["arm"]].add(h)
        i = info[k]
        i["names"].add(r["gene_name"]); i["bt"][r["biotype"]] += 1; i["lens"].append(int(r["length"]))
        i["td"].append(int(r["telomere_distance"]))
    with open(a.out, "w") as fh:
        fh.write("family\tn_arms\tarms\thaplotypes_per_arm\tbiotype\tmedian_length\tmedian_telomere_distance\tmembers\n")
        out = []
        for k, arms in fam.items():
            keep = {x: len(v) for x, v in arms.items() if len(v) >= a.min_frac * len(hap_arm[x])}
            if len(keep) < 2:
                continue
            i = info[k]
            out.append((len(keep), k, ",".join(sorted(keep)), ",".join(f"{x}:{keep[x]}" for x in sorted(keep)),
                        max(i["bt"], key=i["bt"].get), int(st.median(i["lens"])), int(st.median(i["td"])),
                        ",".join(sorted(i["names"]))[:200]))
        for o in sorted(out, key=lambda o: (-o[0], o[1])):
            fh.write("\t".join(map(str, (o[1], o[0]) + o[2:])) + "\n")  # family first, as in the header
    log(f"{len(out)} families on >= 2 arms (min_frac {a.min_frac}); wrote {a.out}")


if __name__ == "__main__":
    main()
