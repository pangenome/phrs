#!/usr/bin/env python3
"""Per-flank verdict from the long-read check of the candidate joins (join_read_check.py --merge output: one row per join
and technology; features closer than the 20 kb side windows are one join, and N runs reach their full contig length).
Per join, over ONT and HiFi (tel/far = reads on the 20 kb on the telomere side / far side of the join):
  spanned        >= 2 reads span it in a technology (SPANNED or SPANNED_LOW_MAPQ)
  end_at_join    no read spans it and reads stop at it with a >= 50% telomeric clip, >= 3 in one technology or >= 2 in
                 both: a chromosome ends at the join
  no_reads_tel   no read spans it, tel <= 2 and far >= 10 in a technology: no read supports the sequence between the join
                 and the telomere next to the join
  no_reads_far   the same with the sides swapped
  not_spanned    none of the above
Per flank: correct (every join spanned); wrong (a join end_at_join or no_reads_tel >= --min-piece bp from the telomere:
the terminal piece is not joined to the rest in any read); wrong_tiny (the same, piece < --min-piece); undecided."""
import argparse, csv, collections

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--checks", required=True)
ap.add_argument("--table", required=True, help="table with flank, hprc_t2t_level, hprc_t2t_gaps, paper_phr_bp")
ap.add_argument("--min-piece", type=int, default=10000)
ap.add_argument("--out", required=True)
a = ap.parse_args()
chk = collections.defaultdict(dict)
for r in csv.DictReader(open(a.checks), delimiter="\t"):
    chk[r["flank"], int(r["start"])][r["tech"]] = r
old = {r["flank"]: r for r in csv.DictReader(open(a.table), delimiter="\t")}
per = collections.defaultdict(list)
for (fl, st), t in chk.items():
    p = fl.endswith("_parm"); ev = []; sp = tstop = False; side = None
    for tech, r in sorted(t.items()):
        L, R = int(r["any_left"]), int(r["any_right"]); tel, far = (L, R) if p else (R, L)
        ev.append(f"{tech}:{r['verdict']}:tel{tel}/far{far}:span{r['span_any_mapq']}:stopTel{r['stop_telomeric']}")
        sp |= int(r["span_any_mapq"]) >= 2
        if tel <= 2 and far >= 10:
            side = "no_reads_tel"
        elif far <= 2 and tel >= 10 and side is None:
            side = "no_reads_far"
    st_tel = [int(r["stop_telomeric"]) for r in t.values() if int(r["span_any_mapq"]) == 0]
    tstop = len(st_tel) == len(t) and (max(st_tel) >= 3 or (len(st_tel) == 2 and min(st_tel) >= 2))
    cls = "spanned" if sp else "end_at_join" if tstop else side or "not_spanned"
    r0 = next(iter(t.values()))
    per[fl].append((int(r0["dist_from_telomere"]), r0["kind"], int(r0["end"]) - int(r0["start"]) + 1, cls, ",".join(ev)))
missing = set(old) - set(per)
if missing:
    raise SystemExit(f"FATAL: {len(missing)} flanks without a read check, e.g. {sorted(missing)[:3]}")
cnt = collections.Counter()
with open(a.out, "w") as fh:
    fh.write("flank\tverdict\tfirst_wrong_join_kb_from_telomere\tjoins(kind@kb:len:class)\thprc_t2t_level\thprc_t2t_gaps\tpaper_phr_bp\tevidence\n")
    for fl in sorted(per):
        v = sorted(per[fl]); bad = [x for x in v if x[3] in ("end_at_join", "no_reads_tel")]
        verdict = ("correct" if all(x[3] == "spanned" for x in v) else "wrong" if any(x[0] >= a.min_piece for x in bad)
                   else "wrong_tiny" if bad else "undecided")
        cnt[verdict] += 1; o = old[fl]
        fh.write("\t".join([fl, verdict, f"{min(x[0] for x in bad) / 1000:g}" if bad else ".",
                            ";".join(f"{k}@{d / 1000:g}:{n}:{c}" for d, k, n, c, _ in v), o["hprc_t2t_level"],
                            o["hprc_t2t_gaps"], o["paper_phr_bp"], " | ".join(e for *_, e in v)]) + "\n")
print(f"{sum(cnt.values())} flanks: {dict(cnt)}; joins: {dict(collections.Counter(x[3] for v in per.values() for x in v))}")
