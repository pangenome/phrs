#!/usr/bin/env python3
"""Inheritance test for a haplotype-specific subtelomeric segment with canonical k-mers (R1-M1 test 3 follow-up).

specific: canonical k-mers of a query interval that occur in none of the exclusion sequences
          (the other ends of the same person and every other window of the same arms in the cohort).
count:    occurrences of those k-mers in reads (BAM/CRAM/FASTA/FASTQ, both strands via canonical form),
          per k-mer and in total, with the read count scanned.
Canonical k-mers are required: a forward-only count misses reverse-complement copies (project rule).
"""
import argparse, gzip, sys, time

COMP = str.maketrans("ACGT", "TGCA")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def canon_kmers(seq, k):
    seq = seq.upper()
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - k + 1):
        f = seq[i:i + k]
        if "N" in f:
            continue
        r = rc[n - i - k:n - i]
        yield f if f < r else r


def read_fasta(path, names=None):
    """yield (name, seq) from a plain or gzipped FASTA; names = optional set to keep"""
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (names is None or name in names):
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (names is None or name in names):
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("specific")
    s.add_argument("--query-fasta", required=True, help="FASTA with the query interval(s), pre-extracted")
    s.add_argument("--exclude-fasta", nargs="+", required=True, help="FASTA files, plain or gzipped")
    s.add_argument("--exclude-names", help="file with sequence names to use from --exclude-fasta (default all)")
    s.add_argument("--k", type=int, default=31)
    s.add_argument("--out", required=True)
    c = sub.add_parser("count")
    c.add_argument("--kmers", required=True)
    c.add_argument("--reads", nargs="+", required=True)
    c.add_argument("--reference", help="FASTA for CRAM decoding")
    c.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "specific":
        qs = list(read_fasta(a.query_fasta))
        q = "N".join(v for _, v in qs)
        qk = set(canon_kmers(q, a.k))
        log(f"query {[n for n, _ in qs]}: {len(q)} bp, {len(qk)} distinct canonical {a.k}-mers")
        names = set(l.strip() for l in open(a.exclude_names)) if a.exclude_names else None
        hit, nseq, nbp = set(), 0, 0
        for path in a.exclude_fasta:
            for name, seq in read_fasta(path, names):
                nseq += 1
                nbp += len(seq)
                hit.update(x for x in canon_kmers(seq, a.k) if x in qk)
        if names is not None and nseq != len(names):
            raise SystemExit(f"FATAL: {len(names)} exclusion names requested, {nseq} found")
        spec = sorted(qk - hit)
        log(f"exclusion: {nseq} sequences, {nbp} bp; query k-mers seen there {len(hit)}; specific {len(spec)}")
        if not spec:
            raise SystemExit("FATAL: no specific k-mers")
        # positions of specific k-mers along the query, for the record
        pos = {}
        for i, x in enumerate(canon_kmers(q, a.k)):
            pos.setdefault(x, i)
        with open(a.out, "w") as fh:
            fh.write("kmer\tquery_offset\n")
            for x in sorted(spec, key=lambda y: pos.get(y, -1)):
                fh.write(f"{x}\t{pos.get(x, 'NA')}\n")
        log(f"wrote {a.out}")
    else:
        import pysam
        km = [l.split("\t")[0] for l in open(a.kmers).read().splitlines()[1:]]
        ks = set(km)
        k = len(km[0])
        with open(a.out, "w") as fh:
            fh.write("reads\tn_reads\tn_specific_kmers\tkmers_seen\tkmer_occurrences\n")
            for path in a.reads:
                cnt = dict.fromkeys(ks, 0)
                n = 0
                af = pysam.AlignmentFile(path, reference_filename=a.reference) if path.endswith((".bam", ".cram")) else None
                it = (r.query_sequence for r in af.fetch(until_eof=True) if r.query_sequence and not r.is_secondary
                      and not r.is_supplementary) if af else (r.sequence for r in pysam.FastxFile(path))
                for seq in it:
                    n += 1
                    for x in canon_kmers(seq, k):
                        if x in ks:
                            cnt[x] += 1
                seen = sum(1 for v in cnt.values() if v)
                occ = sum(cnt.values())
                log(f"{path}: {n} reads, {seen} of {len(ks)} specific k-mers seen, {occ} occurrences")
                fh.write(f"{path}\t{n}\t{len(ks)}\t{seen}\t{occ}\n")


if __name__ == "__main__":
    main()
