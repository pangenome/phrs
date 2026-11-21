#!/usr/bin/env python3
"""Score vote rules on simulated terminal errors.

usage: g1_terminal_rule_check.py WASHU LABEL:GRID_WINDOW_DIR[:CHECK_TSV]... (CHECK_TSV in GRID_WINDOW_DIR, default check.tsv)
For events named *term_* (within 20 kb of a contig end): best ex:5000 candidate under 'more' (HiFi and ONT more)
and 'ont_only' (ONT more); right = NM, length difference <= 50 and matches >= original - 50 (as g1_grid_score.py).
"""
import glob, os, sys
D = sys.argv[1]; os.chdir(D); sys.path.insert(0, "scripts_tmp")
from g1_grid_score_v2 import votes, passes
for spec in sys.argv[2:]:
    label, gd, *cf = spec.split(":")
    chk = {}
    for l in list(open(f"{gd}/{cf[0] if cf else 'check.tsv'}"))[1:]:
        t = l.rstrip("\n").split("\t"); chk[(t[0], t[1])] = (int(t[6]), int(t[7]), t[8], int(t[9]))
    for sf in sorted(glob.glob(f"{gd}/*term_*/summary.tsv")):
        ev = os.path.dirname(sf); name = os.path.basename(ev)
        rows = [l.rstrip("\n").split("\t") for l in open(sf)]; h0 = rows[0]; rows = [dict(zip(h0, r)) for r in rows[1:]]
        v = votes(ev); res = {}
        for rule in ("more", "ont_only"):
            best = None
            for r in rows:
                if r["rule"] != "ex:5000" or r["id"] == "NA": continue
                h, o = v["hifi"][r["id"]], v["ont"][r["id"]]
                ok = passes("more", *o) if rule == "ont_only" else passes("more", *h) and passes("more", *o)
                if ok and (best is None or h[1] + o[1] > best[1]): best = (r["id"], h[1] + o[1], h, o)
            if best is None: res[rule] = "missed"; continue
            ol, mt, nm, ld = chk[(name, best[0])]
            good = nm != "NA" and int(nm) <= 50 and abs(ld) <= 50 and mt >= ol - 50
            res[rule] = "%s %s (NM %s, len diff %d, HiFi %d/%d, ONT %d/%d)" % (best[0], "right" if good else "wrong", nm, ld, best[2][0], best[2][1], best[3][0], best[3][1])
        print(label, name, "| more:", res["more"], "| ont_only:", res["ont_only"])
