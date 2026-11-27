#!/usr/bin/env python3
"""Per-event Flagger status of site-level conversion events, and the counts after removing flagged events.
usage: sites_flagger_events.py FLAGGER_DIR EVENTS_TSV OUT_TSV
  FLAGGER_DIR: HPRC r2 HMM-Flagger v1.2.1 BEDs (<sample>.bed, as conversion_flagger_overlap.py). EVENTS_TSV: events
  (conversion_sites_aggregate.py columns, e.g. sites_cluster_null.py events_significant). An event is FLAGGED when >= 1 bp of
  its span overlaps a flag (same test as conversion_flagger_overlap.py); events of samples without a BED are 'no_bed' and dropped.
out: OUT_TSV = EVENTS_TSV + column flagger (flagged, clean, no_bed); counts of clean events on stderr.
"""
import bisect, collections, csv, glob, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conversion_flagger_overlap import win
FD, EV, OUT = sys.argv[1:4]
flags = collections.defaultdict(list); beds = sorted(glob.glob(f"{FD}/*.bed"))
for f in beds:
    for l in open(f):
        if not l.startswith("track"): c, s, e = l.split("\t")[:3]; flags[c].append((int(s), int(e)))
for v in flags.values(): v.sort()
samples = {os.path.basename(f)[:-4] for f in beds}
def flagged(c, s, e):
    v = flags.get(c, []); i = bisect.bisect_left(v, (e, -1))
    return any(x < e and y > s for x, y in v[max(0, i - 50):i])
rows = list(csv.DictReader(open(EV), delimiter="\t")); n = collections.Counter()
for r in rows:
    c, off, _ = win(r["target"])
    r["flagger"] = "no_bed" if c.split("#")[0] not in samples else "flagged" if flagged(c, off + int(r["start"]), off + int(r["end"])) else "clean"
    n[r["flagger"]] += 1
with open(OUT, "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
cl = [r for r in rows if r["flagger"] == "clean"]
print(f"{len(beds)} BEDs; events {len(rows)}: {dict(n)}", file=sys.stderr)
print(f"clean events {len(cl)}: one individual {sum(r['n_individuals_locus'] == '1' for r in cl)}, shared {sum(r['n_individuals_locus'] not in ('1', 'NA') for r in cl)}, "
      f"outside 2 kb tracts {sum(r['bin_tract_overlap'] == '0' for r in cl)}, embedded {sum(r['cls'] == 'embedded' for r in cl)}, terminal {sum(r['cls'] == 'terminal' for r in cl)}, "
      f"median span {statistics.median(int(r['span']) for r in cl)} bp, median sites {statistics.median(int(r['n_sites']) for r in cl)}", file=sys.stderr)
