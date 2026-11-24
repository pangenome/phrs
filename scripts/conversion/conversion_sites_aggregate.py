#!/usr/bin/env python3
"""Aggregate the site-level conversion scan (conversion_sites.py chunk outputs), R1-M1 test 3b.

Per window: tracts of different donor arms that overlap (first..last site) are merged into one EVENT; its
donors are ranked by the number of sites and the first is the primary donor (a segment copied from arm Y
also carries the bases Y shares with other arms, which gives nested secondary calls). Event classes:
EMBEDDED (sites of the window's own arm on both sides), TERMINAL (open on the telomere side: no own-arm site
between the event and the window end at the telomere, i.e. the end of the chromosome carries the donor
sequence), OPEN (open on the centromeric side of the window, not interpretable). Telomere distance as in
conversion_aggregate.py (p windows start at the telomere, q windows end at it).
LOCUS = events of the same arm and primary donor merged across individuals where their telomere-distance
spans overlap or touch; n_individuals counts samples. C7, C14 and C15 are reported apart (copies there are
as similar as alleles). With --bin-tracts (conversion_aggregate.py <out>.tracts.tsv) we also report how
many events overlap a 2 kb-bin tract of the same window and how many bin tracts contain an event.
OUTPUT: <out>.events.tsv, <out>.loci.tsv, <out>.community.tsv; log on stderr.
"""
import argparse, csv, glob, os, re, statistics as st, sys, time
from collections import defaultdict

NAME = re.compile(r"^([^#]+)#([^#]+)#.*:(\d+)-(\d+)_(chr[0-9XYM]+_(?:parm|qarm))$")


def log(m):
    sys.stderr.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunks", required=True)
    ap.add_argument("--communities", required=True)
    ap.add_argument("--bin-tracts")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    comm = {r["ChromArm"].replace("_p", "_parm").replace("_q", "_qarm"): r["Community"]
            for r in csv.DictReader(open(a.communities), delimiter="\t")}
    tf = sorted(glob.glob(os.path.join(a.chunks, "chunk_*.tracts.tsv")))
    pf = sorted(glob.glob(os.path.join(a.chunks, "chunk_*.power.tsv")))
    if not tf or len(tf) != len(pf):
        raise SystemExit(f"FATAL: {len(tf)} tracts and {len(pf)} power files")
    windows = {r["target"] for f in pf for r in csv.DictReader(open(f), delimiter="\t")}
    by_t = defaultdict(list)
    n = 0
    for f in tf:
        for r in csv.DictReader(open(f), delimiter="\t"):
            n += 1
            by_t[r["target"]].append(r)
    log(f"{len(tf)} chunks, {len(windows)} windows with power rows, {n} tracts in {len(by_t)} windows")
    events = []
    for t, rs in by_t.items():
        m = NAME.match(t)
        if not m:
            raise SystemExit(f"FATAL: cannot parse {t!r}")
        sample, hap, w0, w1, arm = m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)
        L = w1 - w0 + 1
        rs.sort(key=lambda r: int(r["first_site"]))
        cur = []
        for r in rs + [None]:
            if r is not None and cur and int(r["first_site"]) <= max(int(x["last_site"]) for x in cur):
                cur.append(r)
                continue
            if cur:
                s, e = min(int(x["first_site"]) for x in cur), max(int(x["last_site"]) for x in cur) + 1
                don = sorted(cur, key=lambda x: -int(x["n_sites"]))
                p = don[0]
                tel_open = (p["left_a_site"] == "NA") if arm.endswith("parm") else (p["right_a_site"] == "NA")
                cen_open = (p["right_a_site"] == "NA") if arm.endswith("parm") else (p["left_a_site"] == "NA")
                cls = "embedded" if p["embedded"] == "1" else ("open" if cen_open else "terminal")
                d0, d1 = (s, e) if arm.endswith("parm") else (L - e, L - s)
                events.append(dict(target=t, sample=sample, hap=hap, arm=arm, community=comm.get(arm, "none"),
                                   start=s, end=e, span=e - s, d0=d0, d1=d1, cls=cls, primary_donor=p["donor_arm"],
                                   n_sites=int(p["n_sites"]),
                                   donors=",".join(f"{x['donor_arm']}:{x['n_sites']}" for x in don)))
            cur = [r] if r is not None else []
    log(f"{len(events)} events: " + ", ".join(f"{c} {sum(e['cls'] == c for e in events)}"
                                              for c in ("embedded", "terminal", "open")))
    groups = defaultdict(list)
    for i, e in enumerate(events):
        if e["cls"] != "open":
            groups[(e["arm"], e["primary_donor"])].append(i)
    loci = []
    for (arm, don), idx in sorted(groups.items()):
        cur = None
        for i in sorted(idx, key=lambda i: (events[i]["d0"], events[i]["d1"])):
            e = events[i]
            if cur and e["d0"] <= cur["d1"]:
                cur["d1"] = max(cur["d1"], e["d1"]); cur["samples"].add(e["sample"]); cur["ev"].append(i)
            else:
                cur = dict(arm=arm, donor=don, community=e["community"], d0=e["d0"], d1=e["d1"],
                           samples={e["sample"]}, ev=[i])
                loci.append(cur)
    for L_ in loci:
        for i in L_["ev"]:
            events[i]["n_individuals_locus"] = len(L_["samples"])
    if a.bin_tracts:
        bt = defaultdict(list)
        for r in csv.DictReader(open(a.bin_tracts), delimiter="\t"):
            bt[r["target"]].append((int(r["start"]), int(r["end"])))
        hit_e = 0
        for e in events:
            e["bin_tract_overlap"] = int(any(x < e["end"] and y > e["start"] for x, y in bt.get(e["target"], [])))
            hit_e += e["bin_tract_overlap"]
        nb = sum(len(v) for v in bt.values())
        evt = defaultdict(list)
        for e in events:
            evt[e["target"]].append((e["start"], e["end"]))
        hit_b = sum(any(s < y and e > x for s, e in evt.get(t, [])) for t, v in bt.items() for x, y in v)
        log(f"events overlapping a 2 kb-bin tract: {hit_e}/{len(events)}; bin tracts with an event: {hit_b}/{nb}")
    cols = ["target", "sample", "hap", "arm", "community", "start", "end", "span", "d0", "d1", "cls", "primary_donor",
            "n_sites", "donors", "n_individuals_locus"] + (["bin_tract_overlap"] if a.bin_tracts else [])
    with open(a.out + ".events.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in events:
            fh.write("\t".join(str(e.get(c, "NA")) for c in cols) + "\n")
    with open(a.out + ".loci.tsv", "w") as fh:
        fh.write("community\tarm\tprimary_donor\td0\td1\tn_events\tn_individuals\n")
        for x in loci:
            fh.write(f"{x['community']}\t{x['arm']}\t{x['donor']}\t{x['d0']}\t{x['d1']}\t{len(x['ev'])}\t{len(x['samples'])}\n")
    nwin = defaultdict(int)
    for t in windows:
        nwin[comm.get(NAME.match(t).group(5), "none")] += 1
    with open(a.out + ".community.tsv", "w") as fh:
        fh.write("community\twindows\tembedded\tterminal\topen\tloci\tloci_one_individual\tmedian_span\tmedian_sites\n")
        for c in sorted(nwin, key=lambda c: (c == "none", c)):
            ev = [e for e in events if e["community"] == c]
            lc = [x for x in loci if x["community"] == c]
            kept = [e for e in ev if e["cls"] != "open"]
            fh.write(f"{c}\t{nwin[c]}\t" + "\t".join(str(sum(e['cls'] == k for e in ev)) for k in ("embedded", "terminal", "open"))
                     + f"\t{len(lc)}\t{sum(len(x['samples']) == 1 for x in lc)}\t"
                     + (f"{int(st.median(e['span'] for e in kept))}\t{int(st.median(e['n_sites'] for e in kept))}" if kept else "NA\tNA")
                     + "\n")
    log(f"{len(loci)} loci, {sum(len(x['samples']) == 1 for x in loci)} in one individual; wrote {a.out}.*")


if __name__ == "__main__":
    main()
