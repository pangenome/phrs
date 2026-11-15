#!/bin/bash
# Corrected DAD v1.1, own reads only
# run on octopus01: bash g1_build_dad_v11.sh ; out $D/g1/corrected_v11/norm/PAN011#{1,2}.fa.gz + corrections.tsv
# PAN011#2: chr9 end = patch polished with DAD HiFi (patch_v11_hifi/norm), chr7 1,076,488-1,120,704 = hifiasm --ul contig
# (DAD HiFi + DAD ONT UL, g1/local_asm_ul/chr7_h2). PAN011#1 unchanged (copied).
set -euo pipefail
D=/moosefs/guarracino/HPRCv2/PHR_III/pedigrees/washu_v131; cd $D; O=g1/corrected_v11; mkdir -p $O/norm
S=/gnu/store/wzsmrd0x9qb0bdydrg479fq7z5p2v3cc-samtools-1.19/bin/samtools; BGZ=$(ls -d /gnu/store/*htslib-1.19*/bin/bgzip | head -1)
python3 scripts_tmp/g1_apply_patch.py patch_v11_hifi/norm/PAN011#2.fa.gz $O/PAN011#2.fa 'PAN011#2#chr7' 1076488 1120704 \
  g1/local_asm_ul/chr7_h2/contigs.fa loc.bp.hap1.p_ctg.h1tg000001l 228314 272617 + 2> $O/chr7.log
cat $O/chr7.log
$BGZ -@ 4 -c $O/PAN011#2.fa > $O/norm/PAN011#2.fa.gz && $S faidx $O/norm/PAN011#2.fa.gz
cp norm_v11/PAN011#1.fa.gz norm_v11/PAN011#1.fa.gz.fai norm_v11/PAN011#1.fa.gz.gzi $O/norm/
printf "version\tsample\tcontig\tstart\tend\told_bp\tnew_bp\tsource\tevidence\n" > $O/corrections.tsv
printf "v11\tPAN011\tPAN011#2#chr9\t133673491\t133744561\t71071\t69705\tDAD ONT consensus, polished with DAD HiFi\tONT 57 -> 0 wrong, HiFi 21 -> 0, merqury no new k-mers; = DAD-only hifiasm contig (NM 25)\n" >> $O/corrections.tsv
printf "v11\tPAN011\tPAN011#2#chr7\t1076489\t1120704\t44216\t44303\thifiasm --ul (DAD HiFi + DAD ONT UL)\tHiFi depth 2 at the error; contig anchors 206,467 and 201,735 bp at 0.99996 and 0.99990\n" >> $O/corrections.tsv
grep -P "chr7\t|chr9\t" $O/norm/PAN011#2.fa.gz.fai | cut -f1,2; cat $O/corrections.tsv
