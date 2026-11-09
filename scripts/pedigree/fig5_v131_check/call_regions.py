#!/usr/bin/env python3
"""BED of Fig5 call regions on the child diploid, for one version.

Calls come from the ribbon runs of both versions (v1.1 and v1.3.1). Positions
are carried across versions as distance from the nearest chromosome end, then
padded by FLANK. Per call, the region is added on both haplotypes of the query
chromosome and on the partner chromosome end of both haplotypes.

With mode d (donor), the call is placed on the donor side instead: the
target chromosome at donor_start-donor_end, with source lengths taken from the
donor haplotype fai (h1 or h2), and the query chromosome as partner.

usage: call_regions.py OUT.bed TARGET_FAI_H1 TARGET_FAI_H2 SPEC [SPEC ...]
SPEC: q:SRC_FAI:runs.tsv  or  d:SRC_FAI_H1,SRC_FAI_H2:runs.tsv
"""
import csv, sys

FLANK = 1_000_000
ARM = 3_000_000


def load_fai(path):
    return {l.split('\t')[0].split('#')[-1]: int(l.split('\t')[1]) for l in open(path)}


def main():
    out, fai1, fai2, specs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    tgt = {1: load_fai(fai1), 2: load_fai(fai2)}
    hap_name = {h: open(f).readline().split('\t')[0].rsplit('#', 1)[0] for h, f in ((1, fai1), (2, fai2))}
    iv, n_calls = [], 0
    for spec in specs:
        mode, src_fai, runs = spec.split(':')
        srcs = [load_fai(f) for f in src_fai.split(',')]
        for r in csv.DictReader(open(runs), delimiter='\t'):
            n_calls += 1
            if mode == 'q':
                qc, s, e, pc, src = r['query_chrom'], int(r['query_start']), int(r['query_end']), r['target_chrom'], srcs[0]
                d = int(r['donor_start'])
            else:
                qc, s, e, pc = r['target_chrom'], int(r['donor_start']), int(r['donor_end']), r['query_chrom']
                src = srcs[int(r['donor_haplotype'][1]) - 1]
                d = int(r['query_start'])
            qlen = src[qc]
            near_p = s < qlen - e
            for h in (1, 2):
                L = tgt[h].get(qc)
                if L is None:
                    sys.stderr.write(f"WARN {hap_name[h]} lacks {qc}\n"); continue
                if near_p:
                    a, b = s - FLANK, e + FLANK
                else:
                    a, b = L - (qlen - s) - FLANK, L - (qlen - e) + FLANK
                iv.append((f"{hap_name[h]}#{qc}", max(0, a), min(L, b)))
                P = tgt[h].get(pc)
                if P is None:
                    sys.stderr.write(f"WARN {hap_name[h]} lacks partner {pc}\n"); continue
                iv.append((f"{hap_name[h]}#{pc}", 0, min(P, ARM)) if d < P // 2 else (f"{hap_name[h]}#{pc}", max(0, P - ARM), P))
    iv.sort()
    merged = []
    for c, a, b in iv:
        if merged and merged[-1][0] == c and a <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], b)
        else:
            merged.append([c, a, b])
    with open(out, 'w') as o:
        for c, a, b in merged:
            o.write(f"{c}\t{a}\t{b}\n")
    bp = sum(b - a for _, a, b in merged)
    sys.stderr.write(f"{out}: calls={n_calls} intervals={len(merged)} bp={bp}\n")


if __name__ == '__main__':
    main()
