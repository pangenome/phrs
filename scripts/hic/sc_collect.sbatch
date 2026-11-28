#!/bin/bash
#SBATCH -c 2
#SBATCH -t 04:00:00
# Collect, into one small directory, every source file behind the numbers of
# paper_prep/manuscript_revision/similarity_contact_evidence.md that is not yet in the repository, and compute the numbers
# that were first obtained with one-off commands (so that scripts/hic/similarity_contact_numbers.py can print all of
# them from files). Copies: the paper's and every rerun's global_test.tsv (Mantel and within/between), error-rate tables,
# the matched-duplication pairs tables of CHM13 Hi-C (all reads and MAPQ >= 2, 50 kb), the CHM13 duplication-graph
# similarity, the CiFi rerun log (filtered) and its pixel comparison line. Extracts: the evaluate and simulate lines of
# every simulation log, the pair counts of the MAPQ filters. Computes: per CHM13 PHR locus of the matched set, bins with a
# balancing weight and raw cis contacts per million at 50 kb in the paper's Hi-C matrix and in the MAPQ >= 2 matrix;
# read lengths of the sampled real reads (profiles, sampled along whole files); Jaccard distribution of the CHM13
# duplication graph and the chr1 p PHR partners.
# usage: sbatch sc_collect.sbatch OUT_DIR   (new directory)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
O=$1; [ -e $O ] && { echo "ERROR $O exists"; exit 1; }; mkdir -p $O/global_test $O/errrate
B=/moosefs/guarracino/HPRCv2/PHR_III; Q=$B/simhic; M=$B/matched_sd; A=$B/analysis/human/community_based; HC=$B/HiC
PY=/moosefs/guarracino/condatools/porec/1.0-rebuilt/bin/python3
log () { echo "[$(date -Is)] $*"; }
g () { [ -e "$2" ] || { echo "ERROR missing $2"; exit 1; }; cp "$2" $O/global_test/$1.global_test.tsv; }
for r in 50000 20000; do
  for s in hg002 chm13 hg002_porec hg002_cifi; do g paper_${s}_$r $A/${r}bp/${s}_global_test.tsv; done
  for G in HG002 CHM13; do for t in 1 2; do g mapq_${G}_mapq${t}_$r $(ls $Q/mapq/stats_${G}_mapq${t}_$r/*_global_test.tsv); done; done
  g null_hic_hg002_clean_$r $Q/null_pk40/stats_$r/hg002_global_test.tsv
  g null_hic_chm13_clean_$r $Q/chm13_null_rl150/stats_$r/chm13_global_test.tsv
  g null_hic_hg002_realistic_$r $Q/null_wprof_pk40/stats_$r/hg002_global_test.tsv
  g null_hic_chm13_realistic_$r $Q/chm13_null_wprof_pk40/stats_$r/chm13_global_test.tsv
  for mq in 1 0; do
    g null_porec_clean_mq${mq}_$r $Q/porec_null_rk20/stats_mq${mq}_$r/hg002_global_test.tsv
    g null_cifi_clean_mq${mq}_$r $Q/cifi_null_rk20/stats_mq${mq}_$r/hg002_global_test.tsv
    g null_porec_realistic_mq${mq}_$r $Q/porec_null_wprof/stats_mq${mq}_$r/hg002_porec_global_test.tsv
    g rerun_cifi_mq${mq}_$r $Q/rerun_mapq/cifi/stats_mq${mq}_$r/hg002_cifi_global_test.tsv
  done
done
for d in $Q/errrate* ; do if [ -e $d/rates.tsv ]; then cp $d/rates.tsv $O/errrate/$(basename $d).rates.tsv; else log "no rates.tsv in $d (job cancelled before writing)"; fi; done
log "copied $(ls $O/global_test | wc -l) global tests, $(ls $O/errrate | wc -l) error-rate tables"
for f in chm13_hic_50000 chm13_hic_mapq2_50000; do cp $M/contacts/$f.pairs.tsv $O/sd_$f.pairs.tsv; done
cp $(ls $M/pggb_chm13_cen/*.smooth.final.similarity.tsv.gz) $O/sd_graph_chm13.similarity.tsv.gz
cp $Q/logs/rerun_cifi.2887441.log.filtered $O/rerun_cifi.log.filtered
grep -h 'differing pixels' $Q/logs/rerun_followup_cifi.2887971.log | grep -v '^+' > $O/rerun_cifi.compare.txt
# simulate / evaluate lines of every simulation (logs named in the md sections 6 and 9)
for L in simhic.2882447 simhic13b.2887477 simporec.2882915 simcifi.2882916 wsim_hg002.2889804 wsim_chm13.2889810 wsim_porec.2889814; do
  grep -h -E '^(simulate|evaluate|profiles)|min-mapq [01]:|HiC-Pro done' $Q/logs/$L.log | grep -v '^+' | sed "s/^/$L\t/"; done > $O/sim_log_lines.tsv
grep -h -E '^total [0-9]+; MAPQ' $Q/logs/mapqHG002.2882731.log $Q/logs/mapqCHM13.2882732.log | grep -v '^+' > $O/mapq_filter_counts.txt
log "$(wc -l < $O/sim_log_lines.tsv) simulation log lines; MAPQ filter counts: $(cat $O/mapq_filter_counts.txt | tr '\n' ' ')"
$PY - $O $M $HC $Q <<'PY'
import collections, gzip, statistics, sys
import cooler, numpy as np
O, M, HC, Q = sys.argv[1:5]
# CHM13 PHR loci of the matched set: balanced bins and raw cis contacts per million at 50 kb, all reads vs MAPQ >= 2
loci = set()
for l in open(f"{M}/chm13_msd_cen_c0.5.pairs.tsv"):
    x = l.rstrip("\n").split("\t")
    if x[0] == "PHR": loci |= {(x[2], int(x[3]), int(x[4])), (x[7], int(x[8]), int(x[9]))}
C = {k: cooler.Cooler(f"{p}::resolutions/50000") for k, p in (("all", f"{HC}/CHM13/hicpro_output/cool/CHM13.mcool"), ("mapq2", f"{Q}/mapq/CHM13.mapq2.mcool"))}
with open(f"{O}/chm13_phr_loci_50kb.tsv", "w") as o:
    o.write("locus\tlength\tbins\tweighted_bins_all\tweighted_bins_mapq2\traw_cis_per_million_all\traw_cis_per_million_mapq2\n")
    for c, s, e in sorted(loci):
        if c not in C["all"].chromnames or e - s < 5000: continue
        r = f"{c}:{s}-{e}"; row = [r, e - s]
        b = {k: v.bins().fetch(r) for k, v in C.items()}; row.append(len(b["all"]))
        row += [int(np.isfinite(b[k]["weight"].values).sum()) for k in ("all", "mapq2")]
        row += [round(float(v.matrix(balance=False).fetch(r).sum()) / v.info["sum"] * 1e6, 2) for v in C.values()]
        o.write("\t".join(map(str, row)) + "\n")
# read lengths of the real reads sampled along whole files (make_profiles.v2)
with open(f"{O}/read_lengths.tsv", "w") as o:
    o.write("sample\tmate\treads\tmean\tmin\tmax\n")
    for s in ("hg002", "chm13"):
        L = collections.defaultdict(list)
        with gzip.open(f"{Q}/profiles/{s}_hic_whole.tsv.gz", "rt") as f:
            f.readline()
            for l in f:
                x = l.split("\t", 2); L[x[0]].append(int(x[1]))
        for m, v in sorted(L.items()): o.write(f"{s}\t{m}\t{len(v)}\t{statistics.mean(v):.2f}\t{min(v)}\t{max(v)}\n")
# Jaccard distribution of the CHM13 duplication graph and the chr1 p PHR partners
J = []; ex = {}
with gzip.open(f"{O}/sd_graph_chm13.similarity.tsv.gz", "rt") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t")
        if x[0] == x[1]: continue
        J.append(float(x[5]))
        if x[0] == "L1#0#chr1:2-210445" and x[1] in ("L119#0#chr9:150390337-150617247", "L76#0#chr6:171971594-172126628"): ex[x[1]] = x[2:6]
J.sort(); n = len(J)
with open(f"{O}/sd_graph_jaccard.tsv", "w") as o:
    o.write("measure\tvalue\n")
    for k, v in (("ordered_pairs_of_distinct_loci", n), ("median", J[n // 2]), ("p90", J[int(n * 0.9)]), ("p99", J[int(n * 0.99)]), ("max", J[-1])): o.write(f"{k}\t{v}\n")
    for k, v in sorted(ex.items()): o.write(f"chr1p_vs_{k}\tlength_a {v[0]}, length_b {v[1]}, intersection {v[2]}, jaccard {v[3]}\n")
print("done", file=sys.stderr)
PY
ls -la $O $O/global_test | head -80
log "done"
