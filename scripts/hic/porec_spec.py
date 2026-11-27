#!/usr/bin/env python3
"""Exact package list of a damaged conda env: conda packages from conda-meta/*.json (url, md5) plus those whose JSON
is listed in conda-meta/FILES_DELETED.txt (name-version-build from the file name; url and md5 from
PKGS/<pkg>/info/repodata_record.json of the conda package cache); pip packages from site-packages/*.dist-info with
INSTALLER pip or no INSTALLER (name==version). Writes OUT.explicit.txt (@EXPLICIT url#md5), OUT.pip.txt, and lists
the deleted packages not found in the cache.
usage: porec_spec.py ENV PKGS OUT"""
import glob, json, os, sys
env, pkgs, out = sys.argv[1:4]
recs, missing = {}, []
for f in glob.glob(f"{env}/conda-meta/*.json"):
    d = json.load(open(f)); recs[os.path.basename(f)[:-5]] = (d["url"], d["md5"])
deleted = [os.path.basename(l.strip())[:-5] for l in open(f"{env}/conda-meta/FILES_DELETED.txt") if l.strip().endswith(".json")]
for p in deleted:
    r = f"{pkgs}/{p}/info/repodata_record.json"
    if os.path.exists(r):
        d = json.load(open(r)); recs[p] = (d["url"], d["md5"])
    else: missing.append(p)
with open(out + ".explicit.txt", "w") as o:
    o.write("@EXPLICIT\n")
    for p in sorted(recs): o.write(f"{recs[p][0]}#{recs[p][1]}\n")
conda_names = {p.rsplit("-", 2)[0].lower().replace("_", "-") for p in list(recs) + missing}
pip = []
for d in sorted(glob.glob(f"{env}/lib/python3.12/site-packages/*.dist-info")):
    inst = open(f"{d}/INSTALLER").read().strip() if os.path.exists(f"{d}/INSTALLER") else ""
    name, ver = os.path.basename(d)[:-10].rsplit("-", 1)
    if inst != "conda" and name.lower().replace("_", "-") not in conda_names: pip.append(f"{name}=={ver}  # INSTALLER {inst or 'missing'}")
open(out + ".pip.txt", "w").write("\n".join(pip) + "\n")
print(f"conda packages: {len(recs)} ({len(recs) - len(deleted) + len(missing)} from conda-meta, {len(deleted) - len(missing)} deleted JSON found in the cache); "
      f"deleted JSON not in the cache: {missing}; pip packages: {len(pip)}")
