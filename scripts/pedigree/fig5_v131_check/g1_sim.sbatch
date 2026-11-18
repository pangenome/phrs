#!/bin/bash
#SBATCH -p workers,tux
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 4:00:00
#SBATCH --exclude=octopus05
# Assembly with known errors (G1r)
# usage: sbatch g1_sim.sbatch --washu D --sample SM --version V --hap H --seed N --out DIR
# out (new dir only): DIR/norm/SM#{1,2}.fa.gz (+fai, gzi; HAP modified by g1_sim.py, other copied), DIR/truth.tsv,
#   DIR/locus_tasks.tsv (id, name, error on DIR/norm, read error on norm_V; 10 kb bins as the detection gives)
set -euxo pipefail
trap 'echo "exit status $? at line ${LINENO}"' EXIT
while [ $# -gt 0 ]; do case $1 in --washu) D=$2;; --sample) SM=$2;; --version) V=$2;; --hap) H=$2;; --seed) SD=$2;; --out) O=$2;; *) echo "ERROR unknown argument $1"; exit 2;; esac; shift 2; done
for v in D SM V H SD O; do [ -n "${!v:-}" ] || { echo "ERROR missing argument for $v"; exit 2; }; done
log () { echo "[$(date -Is)] $*"; }
cd $D; [ -s $O/truth.tsv ] && { log "exists $O/truth.tsv, skip"; exit 0; }
[ -e $O ] && { echo "ERROR $O exists without truth.tsv; use a new dir"; exit 1; }
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
OH=$((3 - H)); log "$SM $V hap $H seed $SD -> $O"
python3 scripts_tmp/g1_sim.py $D $SM $V $H $SD $O
$BGZ -@ 4 $O/norm/$SM#$H.fa && $S faidx $O/norm/$SM#$H.fa.gz
cp norm_$V/$SM#$OH.fa.gz norm_$V/$SM#$OH.fa.gz.fai norm_$V/$SM#$OH.fa.gz.gzi $O/norm/
log "lengths changed: $(paste <(cut -f1,2 norm_$V/$SM#$H.fa.gz.fai) <(cut -f1,2 $O/norm/$SM#$H.fa.gz.fai) | awk '$2 != $4' | wc -l) contigs"
awk -F'\t' -v OFS='\t' 'NR>1 {b=10000; ms=int($6/b)*b+1; me=($7>$6 ? int(($7+b-1)/b)*b : ms+b-1); os=int($4/b)*b+1; oe=int(($5+b-1)/b)*b; print NR-1, $1"_"$2, $3":"ms"-"me, $3":"os"-"oe}' $O/truth.tsv > $O/locus_tasks.tsv
log "locus tasks $(wc -l < $O/locus_tasks.tsv) (error and read intervals on 10 kb bins)"
cat $O/truth.tsv
log done
