# Subtelomeric PHR detection pipeline

Command-level record of how the subtelomeric pseudo-homolog regions (PHRs) were
detected, from telomere-bearing contig selection through the pangenome graph.
Copied from the source-of-truth lab notes
(`Dropbox/working/Garrison/hprcv2/code_PHR-and-3D.md`, section `## PHR III`) and
the two custom scripts in this directory.

These commands ran on the moosefs deploy environment; absolute `/moosefs/...`
paths and SLURM `sbatch` submissions are reproduced verbatim. `$dir_base` is the
analysis root (`/moosefs/guarracino/HPRCv2`). Tool versions are in the paper's
"Software versions" Methods subsection (PGGB v0.7.4 bundle, wfmash v0.23.0,
impg 0.4.1, odgi v0.9.2).

Stages:
1. Contig classification and telomere filtering (telomere calling + T2T filtering)
2. Subtelomeric flank extraction and telomere trimming (500 kb flanks)
3. All-vs-all flank alignment (wfmash)
4. Inter-chromosomal region detection (PHR detection)
5. Pangenome graph construction and Jaccard similarity (pggb + odgi)
6. Community detection (Leiden) — next stage, included for completeness

> IMPORTANT (transitivity). PHR detection in stage 4 runs `impg query` **without**
> `-x`/`--transitive`, i.e. a direct all-vs-all projection through the
> all-against-all alignments, **not** a transitive closure. impg's `query`
> transitive flag defaults to off (`impg` `src/main.rs`,
> `#[clap(short = 'x', long, action)] transitive: bool`). Transitive closure was
> used for the genome-wide Fig. 1a coverage scan, a separate whole-genome impg
> index, not for the subtelomeric PHR pipeline here.

---

## Stage 1. Contig classification and telomere filtering

Telomere annotations are precomputed per assembly (`*.telo.tsv` under
`/moosefs/pangenomes/HPRCv2/telomeres`). This stage maps every assembly to CHM13,
builds p-arm / q-arm reference intervals, and keeps only contigs >= 1 Mb that
carry a real telomere on the correct end (`classify_contigs.py`).

```shell
# Map all HPRCv2 assemblies to CHM13 (mapping only, -m)
mkdir -p $dir_base/PHR_III/assembly-vs-chm13
cd $dir_base/PHR_III/assembly-vs-chm13
ls /moosefs/pangenomes/HPRCv2/*.fa.gz | while read fasta; do
    sbatch -p allnodes -c 48 -x tux08 --job-name $(basename $fasta)-vs-chm13 --wrap \
      "hostname; cd /scratch; /moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash \
        /moosefs/pangenomes/HPRCv2/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz $fasta \
        -p 90 -t 48 --quiet -m > /scratch/$(basename $fasta).vs.chm13.map.paf; \
       mv /scratch/$(basename $fasta).vs.chm13.map.paf $dir_base/PHR_III/"
done

# p-arm / q-arm reference intervals from CHM13 centromere coordinates
mkdir -p $dir_base/PHR_III/pq-classification
cd $dir_base/PHR_III/pq-classification
grep CHM13 *tsv -h | cut -f 4-6 | sort -V | sed 's/chr/CHM13#0#chr/g' > $dir_base/data/chm13.centromeres.bed
bedtools complement \
    -i $dir_base/data/chm13.centromeres.bed \
    -g <(cut -f 1,2 /moosefs/pangenomes/HPRCv2/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz.fai) | grep chrM -v > tmp.bed
sed -n 1~2p tmp.bed > p_arms.bed   # odd rows  = p-arms
sed -n 2~2p tmp.bed > q_arms.bed   # even rows = q-arms
rm tmp.bed

# Classify contigs: p-contig / q-contig / pq-contig with real telomere on the correct end
python3 classify_contigs.py \
    --p-arms p_arms.bed \
    --q-arms q_arms.bed \
    --chrom-alias /moosefs/guarracino/HPRCv2/coverage/hprcv2.chromAlias.txt \
    --paf-pattern '/moosefs/guarracino/HPRCv2/PHR_III/assembly-vs-chm13/*.paf' \
    --telomere-dir /moosefs/pangenomes/HPRCv2/telomeres \
    --t2t-dir /moosefs/guarracino/HPRCv2/coverage/t2t-assignments \
    --min-contig-length 1000000 \
    --min-alignment-length 0 \
    --output contig_classifications.tsv 2> classification.log
```

`classify_contigs.py` keeps a contig only if it is >= 1 Mb, maps to its own
chromosome's p- and/or q-arm, and carries the matching telomere from the
`.telo.tsv` annotation (p-contig -> p-telomere, q-contig -> q-telomere, reverse
strand swaps them); contigs with duplicate same-type telomeres (assembly errors)
are dropped.

## Stage 2. Subtelomeric flank extraction and telomere trimming

For every valid telomere, extract the telomere plus 500 kb of adjacent sequence,
then **trim the telomere repeat**, keeping only the 500 kb flank internal to it.
`FLANK=500000`. GRCh38 and `CHM13#0#chrY` are excluded here. Output:
`hprcv2.1Mb.telo_500kb_trimmed.fa.gz`.

```shell
cd $dir_base/PHR_III/pq-classification

TELO_DIR="/moosefs/pangenomes/HPRCv2/telomeres"
FASTA_DIR="/moosefs/pangenomes/HPRCv2"
CHROM_ALIAS="/moosefs/guarracino/HPRCv2/coverage/hprcv2.chromAlias.txt"
CLASSIFICATION_FILE="/moosefs/guarracino/HPRCv2/PHR_III/pq-classification/contig_classifications.tsv"
OUTPUT_DIR="/scratch"
FINAL_DIR="/moosefs/guarracino/HPRCv2/PHR_III"
MIN_LEN=1000000
FLANK=500000
THREADS=48
OUTPUT_NAME="hprcv2.1Mb.telo_500kb_trimmed"

# Valid contigs from stage 1
tail -n +2 "${CLASSIFICATION_FILE}" | cut -f1 > "${OUTPUT_DIR}/valid_contigs.txt"

# Telomere records, excluding GRCh38 and CHM13 chrY, restricted to valid contigs
cat "${TELO_DIR}"/*.telo.tsv | grep -v "^GRCh38" | grep -v "^CHM13#0#chrY" | \
    grep -Ff "${OUTPUT_DIR}/valid_contigs.txt" > "${OUTPUT_DIR}/telomeres.tsv"

# extract.bed = telomere + 500 kb ; trim.bed = the 500 kb flank only (telomere removed)
rm -f "${OUTPUT_DIR}/extract.bed" "${OUTPUT_DIR}/trim.bed"
awk -F'\t' -v flank="$FLANK" -v min_len="$MIN_LEN" -v outdir="$OUTPUT_DIR" '
BEGIN {OFS="\t"}
NR==FNR {chrom[$1]=$2; next}
{
    seq_id = $1; start = $2; end = $3; seq_len = $4
    if (seq_len < min_len) next
    chr_suffix = (seq_id in chrom) ? "_" chrom[seq_id] : ""
    if (start == 0) {                                   # p-arm telomere at the start
        extract_end = (end + flank > seq_len) ? seq_len : end + flank
        name = seq_id ":1-" extract_end chr_suffix "_parm"
        print seq_id, 0, extract_end, name >> outdir"/extract.bed"
        print name, end, extract_end >> outdir"/trim.bed"
    }
    if (end == seq_len) {                               # q-arm telomere at the end
        extract_start = (start - flank < 0) ? 0 : start - flank
        name = seq_id ":" (extract_start + 1) "-" seq_len chr_suffix "_qarm"
        print seq_id, extract_start, seq_len, name >> outdir"/extract.bed"
        print name, 0, start - extract_start >> outdir"/trim.bed"
    }
}' "$CHROM_ALIAS" "${OUTPUT_DIR}/telomeres.tsv"

# contig -> source FASTA map
for fasta in "${FASTA_DIR}"/*.fa.gz; do
    [[ -f "${fasta}.fai" ]] && awk -v f="$fasta" '{print $1"\t"f}' "${fasta}.fai"
done > "${OUTPUT_DIR}/contig_map.tsv"

# Extract telomere+flank
awk -F'\t' 'NR==FNR {fasta[$1]=$2; next} {contig=$1;start=$2;end=$3;name=$4;
    if (contig in fasta) print fasta[contig], contig":"(start+1)"-"end, name}' \
    "${OUTPUT_DIR}/contig_map.tsv" "${OUTPUT_DIR}/extract.bed" > "${OUTPUT_DIR}/extraction_commands.tsv"
rm -f "${OUTPUT_DIR}/extracted.fa"
while read -r fasta_file region name; do
    samtools faidx "$fasta_file" "$region" 2>/dev/null | awk -v n="$name" 'NR==1 {print ">"n; next} {print}'
done < "${OUTPUT_DIR}/extraction_commands.tsv" > "${OUTPUT_DIR}/extracted.fa"
samtools faidx "${OUTPUT_DIR}/extracted.fa"

# Trim the telomere: keep only the flank (trim.bed sub-region of each extracted sequence)
rm -f "${OUTPUT_DIR}/${OUTPUT_NAME}.fa"
while read -r seqname start end; do
    base_id=$(echo "$seqname" | sed 's/:[0-9]*-[0-9]*_.*$//')
    orig_start=$(echo "$seqname" | sed 's/.*:\([0-9]*\)-[0-9]*_.*$/\1/')
    suffix=$(echo "$seqname" | sed 's/.*:[0-9]*-[0-9]*\(_.*\)$/\1/')
    new_start=$((orig_start + start)); new_end=$((orig_start + end - 1))
    new_name="${base_id}:${new_start}-${new_end}${suffix}"
    samtools faidx "${OUTPUT_DIR}/extracted.fa" "${seqname}:$((start + 1))-${end}" | \
        awk -v n="$new_name" 'NR==1 {print ">"n; next} {print}'
done < "${OUTPUT_DIR}/trim.bed" > "${OUTPUT_DIR}/${OUTPUT_NAME}.fa"

bgzip -@ "$THREADS" -l 9 "${OUTPUT_DIR}/${OUTPUT_NAME}.fa"
samtools faidx "${OUTPUT_DIR}/${OUTPUT_NAME}.fa.gz"
mv "${OUTPUT_DIR}/${OUTPUT_NAME}.fa.gz"* "${FINAL_DIR}/"
```

This yields the 18,827 telomere-anchored 500 kb flanks (telomere repeat trimmed).

## Stage 3. All-vs-all flank alignment

wfmash all-against-all at 95% identity, each flank serving in turn as the target
(`--target-prefix`), one SLURM array task per target sequence.

```shell
mkdir -p $dir_base/PHR_III/all-vs-all.1Mb.p95
cd $dir_base/PHR_III/all-vs-all.1Mb.p95
cut -f 1 $dir_base/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai > sequences.txt
n_seqs=$(wc -l < sequences.txt)
export dir_base
sbatch -c 96 -p tux -x tux09 --array=1-${n_seqs}%20 --job-name=wfmash_array <<'EOF'
#!/bin/bash
seq=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sequences.txt)
seq_sanitized=$(echo $seq | tr ':@#-' '____')
cd /scratch
\time -v /moosefs/guarracino/pggb_wfmash023/wfmash/build/bin/wfmash \
    -p 95 --target-prefix "$seq" -t 48 \
    $dir_base/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz --quiet \
  | bgzip -l 9 -@ 48 > /scratch/hprcv2.1Mb.telo_500kb_trimmed.all-vs-${seq_sanitized}.paf.gz
mv /scratch/hprcv2.1Mb.telo_500kb_trimmed.all-vs-${seq_sanitized}.paf.gz $dir_base/PHR_III/all-vs-all.1Mb.p95/
EOF
```

## Stage 4. Inter-chromosomal region detection (PHR detection)

Build an impg index over the all-vs-all PAFs, then walk each flank inward from the
telomere in 5 kb windows. A window "passes" when the **direct** (non-transitive)
`impg query` projection reaches >= 2 distinct chromosomes each with >= 5 alignments
at >= 95% identity, with each projected alignment >= 3 kb. Contiguous passing
windows form a PHR; scanning stops after 4 consecutive failing windows.

```shell
cp -r $dir_base/PHR_III/all-vs-all.1Mb.p95 /scratch/all-vs-all.1Mb.p95
impg_index=/scratch/all-vs-all.1Mb.p95.impg
paf_list="$dir_base/PHR_III/all-vs-all.1Mb.p95.paf.list"
find /scratch/all-vs-all.1Mb.p95 -name '*.paf.gz' > "$paf_list"
impg index --alignment-list "$paf_list" -i $impg_index -t 32

cd $dir_base/PHR_III
python3 find-multichr-regions-incremental.py \
    --paf-list "$paf_list" -i $impg_index \
    --fai $dir_base/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz.fai \
    --window 5000 --step 5000 --max 500000 \
    --min-identity 0.95 --min-output-length 3000 \
    --min-diff-chrs 2 \
    --min-count 5 \
    --min-consecutive-stop 4 \
    --per-window-output $dir_base/PHR_III/all-vs-all.1Mb.p95.id95.len.per_window.tsv \
    --threads 40 --verbose \
    > $dir_base/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv

# Drop the NA18982 hap1 chr18_q scaffolding chimera (chr18 joined to chrX PAR1 across a 100 bp N-gap)
grep -v 'NA18982#1#JBKABS010000018.1.*chr18_qarm' \
    $dir_base/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv > tmp && mv tmp $dir_base/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv

# Cut each flank down to its called PHR interval -> the 15,668 PHR sequences
fasta="$dir_base/PHR_III/hprcv2.1Mb.telo_500kb_trimmed.fa.gz"
regions_tsv="$dir_base/PHR_III/all-vs-all.1Mb.p95.id95.len.tsv"
output="$dir_base/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa"
rm -f "$output"
sed 1d "$regions_tsv" | cut -f 1,4,5 | awk '$2 != "."' | while read -r seq_name int_start int_end; do
    prefix=$(echo "$seq_name" | cut -d ':' -f 1)
    coords_suffix=$(echo "$seq_name" | cut -d ':' -f 2)
    orig_start=$(echo "$coords_suffix" | cut -d '-' -f 1)
    orig_end=$(echo "$coords_suffix" | cut -d '-' -f 2 | cut -d '_' -f 1)
    suffix=$(echo "$coords_suffix" | cut -d '-' -f 2 | sed 's/^[0-9]*//')
    if [[ "$seq_name" == *"_parm"* ]]; then
        new_start=$orig_start; new_end=$((orig_start + int_end))
    else
        new_start=$((orig_start + int_start)); new_end=$orig_end
    fi
    new_name="${prefix}:${new_start}-${new_end}${suffix}"
    samtools faidx "$fasta" "${seq_name}:$((int_start + 1))-${int_end}" | sed "1s/.*/>$new_name/" >> "$output"
done
bgzip -@ 24 -l 9 "$output"; samtools faidx "$output.gz"
```

Result: 15,668 PHRs (83.2% of the 18,827 flanks) on 41 of 48 arms. See
`find-multichr-regions-incremental.py` in this directory for the exact
window/threshold logic.

## Stage 5. Pangenome graph and Jaccard similarity

pggb on the 15,668 PHR sequences at 95% identity, then odgi similarity (Jaccard =
fraction of shared graph nodes between each pair of sequences).

```shell
mkdir -p $dir_base/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95
sbatch -w octopus02 -c 48 -p allnodes --job-name=pggb-subtelo-1Mb --wrap \
  "pggb --resume -i $dir_base/PHR_III/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz \
        -o $dir_base/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95 -p 95 -D /scratch -t 48"

dir_graph=$dir_base/PHR_III/pggb/hprcv2.1Mb.telo_trimmed.p95.id95
graph=$dir_graph/hprcv2.1Mb.telo_trimmed.p95.id95.fa.gz.*.smooth.final.og
name=$(basename $graph .og)
sbatch -c 48 -p allnodes --job-name=odgi-sim-1Mb \
  --wrap "odgi similarity -i $graph -t 48 -P --all | pigz -9 > $dir_graph/\$name.similarity.tsv.gz"
```

## Stage 6. Community detection (next stage)

Arm-level (41 arms -> 15 Leiden communities) and sequence-level partitions from
the Jaccard matrix (`detect_communities.R`, in the lab `scripts/community/`; not
copied here). Decompress the similarity TSV to disk first (never pipe zcat into
fread; it truncates silently).

```shell
Rscript detect_communities.R --similarity $SIM_TSV --organism human --level arm \
    --output-prefix $dir_base/PHR_III/similarity/hprcv2.1Mb.subtelo
```
