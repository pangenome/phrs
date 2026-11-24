#!/bin/bash
# Read support for a random sample of conversion-scan tracts (R1-M1 test 3b control): HPRC release-2 Flagger
# HiFi BAMs (reads aligned to the same diploid assembly) and the assembly sequence, both over HTTPS.
# Usage: tract_read_support.sh TRACTS_TSV HPRC_INDEX_DIR N SEED CLASS OUTDIR
#   TRACTS_TSV      conversion_aggregate.py <out>.tracts.tsv
#   HPRC_INDEX_DIR  assemblies_release2_v1.0.index.csv and flagger_hifi_processing_metadata_v0.1.csv
#   CLASS           singleton (n_other_individuals = 0) or shared; C7, C14 and C15 are excluded
# Output: OUTDIR/<i>.tsv per tract (read_support_region.py) and OUTDIR/sample.tsv (the tracts drawn)
set -euxo pipefail
T=$1; IDX=$2; N=$3; SEED=$4; CLASS=$5; O=$6
D=$(cd "$(dirname "$0")" && pwd)
mkdir -p $O
python3 - "$T" "$IDX" "$N" "$SEED" "$CLASS" > $O/sample.tsv <<'EOF'
import csv, random, re, sys
t, idx, n, seed, cls = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
asm = {(r["sample_id"], r["haplotype"]): r["assembly"] for r in csv.DictReader(open(f"{idx}/assemblies_release2_v1.0.index.csv"))}
bam = {r["sample_id"]: r["bam"] for r in csv.DictReader(open(f"{idx}/flagger_hifi_processing_metadata_v0.1.csv"))}
https = lambda s: s.replace("s3://human-pangenomics/", "https://s3-us-west-2.amazonaws.com/human-pangenomics/")
rows = [r for r in csv.DictReader(open(t), delimiter="\t") if r["community"] not in ("C7", "C14", "C15")
        and (r["n_other_individuals"] == "0") == (cls == "singleton") and r["sample"] in bam
        and (r["sample"], r["hap"]) in asm]
sys.stderr.write(f"{len(rows)} {cls} tracts with BAM and assembly; drawing {n} (seed {seed})\n")
for r in random.Random(seed).sample(rows, n):
    m = re.match(r"^(.+):(\d+)-(\d+)_chr", r["target"])
    off = int(m.group(2)) - 1
    print(r["sample"], r["hap"], r["community"], r["arm"], r["paralog_arm"], m.group(1), off + int(r["start"]),
          off + int(r["end"]), https(bam[r["sample"]]), https(asm[(r["sample"], r["hap"])]), sep="\t")
EOF
i=0
while IFS=$'\t' read -r s h c arm par contig st en B A; do
    i=$((i + 1))
    lo=$((st - 30000 > 0 ? st - 30000 : 0)); hi=$((en + 30000))
    W=$(mktemp -d)
    samtools view -b -o $W/r.bam "$B" "$contig:$((lo + 1))-$hi"
    samtools index $W/r.bam
    samtools faidx "$A" "$contig:$((lo + 1))-$hi" > $W/a.fa
    samtools faidx $W/a.fa
    printf '# %s\t%s\t%s\t%s\tdonor %s\t%s:%s-%s\n' $s $h $c $arm $par "$contig" $st $en > $O/$i.tsv
    python3 $D/read_support_region.py --bam $W/r.bam --fasta $W/a.fa --contig "$contig" --start $st --end $en \
        --pad 1000 --flank 10000 >> $O/$i.tsv
    rm -rf $W
done < $O/sample.tsv
