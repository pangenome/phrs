#!/usr/bin/env python3
"""Empirical error profiles of real reads, for null simulations that look like the real data (R1-M4 / R3-M3).

hic   R1_SAM R2_SAM OUT_TSV_GZ  bowtie2 alignments of real Hi-C R1 and R2 reads (the paper's HiC-Pro options). For each
      primary mapped read, in read orientation: mate, length, base qualities, bases soft-clipped at the 5' and 3' ends,
      the read positions of mismatches (from the MD tag) and the read name. Insertions and deletions are not kept (gap extensions are
      0.05-0.08% per aligned base, against 1.5-2.4% mismatches).
porec SAM OUT_TSV_GZ            minimap2 alignments of digested real Pore-C / CiFi monomers (names readID:...). Per read
      with at least one primary mapped monomer: aligned bases (CIGAR M), substitution, insertion and deletion rates per
      aligned base (mismatches = NM - inserted - deleted bases), and per mapped monomer the fractions of its length
      soft-clipped at its 5' and 3' ends (read orientation).
Logs the reads kept and dropped with the reason and the mean rates."""
import collections, gzip, sys
import pysam
log = lambda *m: print(*m, file=sys.stderr, flush=True)
mode = sys.argv[1]
if mode == "hic":
    drop, n, mm, al, cl, bases = collections.Counter(), 0, 0, 0, 0, 0
    with gzip.open(sys.argv[4], "wt") as o:
        o.write("mate\tlength\tqual\tclip5\tclip3\tmismatches\tread\n")
        for mate, p in ((1, sys.argv[2]), (2, sys.argv[3])):
            for r in pysam.AlignmentFile(p):
                if r.is_unmapped: drop["unmapped"] += 1; continue
                if r.is_secondary or r.is_supplementary: drop["secondary or supplementary"] += 1; continue
                L = r.query_length; c = r.cigartuples
                left = c[0][1] if c[0][0] == 4 else 0; right = c[-1][1] if c[-1][0] == 4 else 0
                q = "".join(chr(x + 33) for x in r.query_qualities)
                pos = [qp for qp, rp, rb in r.get_aligned_pairs(matches_only=True, with_seq=True) if rb.islower()]
                if r.is_reverse: q, left, right, pos = q[::-1], right, left, [L - 1 - x for x in pos]
                o.write(f"{mate}\t{L}\t{q}\t{left}\t{right}\t{','.join(map(str, sorted(pos))) or '.'}\t{r.query_name}\n")
                n += 1; mm += len(pos); al += L - left - right; cl += left + right; bases += L
    log(f"hic profiles: {n} reads kept; dropped {dict(drop)}; mismatches per aligned base {mm / al:.5f}; soft-clipped per read base {cl / bases:.5f}")
elif mode == "porec":
    R = collections.defaultdict(lambda: [0, 0, 0, 0, []]); drop = collections.Counter()
    for r in pysam.AlignmentFile(sys.argv[2]):
        if r.is_unmapped: drop["unmapped monomer"] += 1; continue
        if r.is_secondary or r.is_supplementary: drop["secondary or supplementary"] += 1; continue
        c = r.cigartuples; L = r.query_length
        m = sum(x for op, x in c if op == 0); i = sum(x for op, x in c if op == 1); d = sum(x for op, x in c if op == 2)
        left = c[0][1] if c[0][0] == 4 else 0; right = c[-1][1] if c[-1][0] == 4 else 0
        if r.is_reverse: left, right = right, left
        x = R[r.query_name.split(":")[0]]; x[0] += m; x[1] += r.get_tag("NM") - i - d; x[2] += i; x[3] += d; x[4].append(f"{left / L:.4f}/{right / L:.4f}")
    with gzip.open(sys.argv[3], "wt") as o:
        o.write("read\taligned\tsub\tins\tdel\tmonomer_clips\n")
        for k, (m, s, i, d, cl) in R.items():
            o.write(f"{k}\t{m}\t{s / m:.6f}\t{i / m:.6f}\t{d / m:.6f}\t{','.join(cl)}\n")
    M = sum(v[0] for v in R.values())
    log(f"porec profiles: {len(R)} reads, {sum(len(v[4]) for v in R.values())} mapped monomers; dropped {dict(drop)};"
        f" per aligned base sub {sum(v[1] for v in R.values()) / M:.5f}, ins {sum(v[2] for v in R.values()) / M:.5f}, del {sum(v[3] for v in R.values()) / M:.5f}")
else:
    sys.exit(f"FATAL mode {mode}")
