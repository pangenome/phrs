#!/usr/bin/env python3
"""Read support for an assembled interval (candidate exchange or conversion tract), per BAM of reads aligned
to the same diploid assembly (e.g. the HPRC release-2 Flagger HiFi and ONT BAMs).

For each BAM: (1) primary reads that span the whole interval and both flanks of --pad bp, with MAPQ and
the number of X (=/X CIGAR) and indel events inside the interval; (2) pileup over the interval and each
flank (MAPQ >= --min-mapq): median depth and the number of columns where > 20% of >= 3 read bases differ
from the assembly. Needs the assembly sequence of the contig (FASTA, same contig name as the BAM).
"""
import argparse, sys, time
import pysam


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", nargs="+", required=True)
    ap.add_argument("--fasta", required=True, help="assembly FASTA holding --contig")
    ap.add_argument("--contig", required=True)
    ap.add_argument("--start", type=int, required=True, help="0-based")
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--pad", type=int, default=3000)
    ap.add_argument("--flank", type=int, default=10000)
    ap.add_argument("--min-mapq", type=int, default=20)
    a = ap.parse_args()
    fa = pysam.FastaFile(a.fasta)
    off = 0
    name = a.contig
    if name not in fa.references:  # FASTA extracted as contig:1-N by samtools faidx
        hits = [r for r in fa.references if r.startswith(a.contig + ":")]
        if len(hits) != 1:
            raise SystemExit(f"FATAL: {a.contig} not in {a.fasta}")
        name = hits[0]
        off = int(name.rsplit(":", 1)[1].split("-")[0]) - 1
    seq = fa.fetch(name)
    s0, e0 = a.start - a.pad, a.end + a.pad
    log(f"interval {a.contig}:{a.start}-{a.end}, span test {s0}-{e0}, flanks {a.flank} bp, MAPQ >= {a.min_mapq}")
    for path in a.bam:
        bam = pysam.AlignmentFile(path)
        span = []
        for r in bam.fetch(a.contig, a.start, a.end):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.reference_start <= s0 and r.reference_end >= e0:
                x = ev = 0
                pos = r.reference_start
                for op, n in r.cigartuples:
                    if op in (7, 8):
                        if op == 8:
                            x += max(0, min(pos + n, a.end) - max(pos, a.start))
                        pos += n
                    elif op in (1, 2):
                        ev += int(a.start <= pos < a.end)
                        pos += n if op == 2 else 0
                    elif op == 0:
                        raise SystemExit("FATAL: M in CIGAR; need =/X")
                span.append((r.mapping_quality, r.reference_start, r.reference_end, x, ev))
        q = [s for s in span if s[0] >= a.min_mapq]
        print(f"{path}\tspanning_reads\t{len(span)}\tMAPQ>={a.min_mapq}\t{len(q)}\tX_in_interval\t"
              + ",".join(str(s[3]) for s in sorted(q, key=lambda s: s[3])))
        for lab, lo, hi in (("left", a.start - a.flank, a.start), ("interval", a.start, a.end),
                            ("right", a.end, a.end + a.flank)):
            dep, bad, n = [], 0, 0
            for col in bam.pileup(a.contig, max(0, lo), hi, truncate=True, min_mapping_quality=a.min_mapq,
                                  min_base_quality=0, stepper="nofilter"):
                bases = [p.alignment.query_sequence[p.query_position] for p in col.pileups
                         if not p.is_del and not p.is_refskip and p.query_position is not None]
                n += 1
                dep.append(len(col.pileups))
                ref = seq[col.reference_pos - off].upper()
                if len(bases) >= 3 and sum(b != ref for b in bases) / len(bases) > 0.2:
                    bad += 1
            dep.sort()
            print(f"{path}\tpileup_{lab}\t{lo}-{hi}\tmedian_depth\t{dep[len(dep) // 2] if dep else 0}\t"
                  f"discordant_columns\t{bad}\tcolumns\t{n}")


if __name__ == "__main__":
    main()
