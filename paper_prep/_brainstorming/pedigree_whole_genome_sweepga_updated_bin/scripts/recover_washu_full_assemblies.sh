#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${RECOVERY_DIR:-/moosefs/erikg/phrs/recovery/fig5-whole-genome-joint-parent-sweepga}"
BASE_URL="${WASHU_PUBLIC_BASE_URL:-https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1}"
THREADS="${BGZIP_THREADS:-4}"
SAMPLES=("$@")
if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    SAMPLES=(PAN010 PAN011 PAN027 PAN028)
fi

mkdir -p "$OUT_DIR"
manifest="$OUT_DIR/recovery_manifest.tsv"
if [[ ! -s "$manifest" ]]; then
    printf "sample\tsource_url\toutput_fasta\toutput_fai\toutput_gzi\theader_rule\tstatus\n" > "$manifest"
fi

for sample in "${SAMPLES[@]}"; do
    url="$BASE_URL/assembly.v1.1.${sample}.diploid.fa"
    out="$OUT_DIR/${sample}.fa.gz"
    tmp="$out.tmp"
    rm -f "$tmp"
    echo "recovering sample=$sample url=$url out=$out"
    curl -L "$url" |
        awk -v sample="$sample" '
            BEGIN { FS="[ .]" }
            /^>/ && $1 == ">" sample {
                chrom=$2
                hap=$3
                sub(/^haplotype/, "", hap)
                if (hap == "maternal") {
                    hap="1"
                } else if (hap == "paternal") {
                    hap="2"
                }
                print ">" sample "#" hap "#" chrom
                next
            }
            { print }
        ' |
        bgzip -@ "$THREADS" -c > "$tmp"
    mv "$tmp" "$out"
    samtools faidx "$out"
    gzip -t "$out"
    test -s "$out.fai"
    test -s "$out.gzi"
    grep -q "^${sample}#1#chr1" "$out.fai"
    grep -q "^${sample}#2#chr1" "$out.fai"
    awk -v sample="$sample" 'BEGIN{FS=OFS="\t"} $1 != sample {print}' "$manifest" > "$manifest.tmp"
    printf "%s\t%s\t%s\t%s\t%s\t%s\tOK\n" \
        "$sample" "$url" "$out" "$out.fai" "$out.gzi" \
        "PAN sample public headers converted from SAMPLE.chrN.haplotypeH or SAMPLE.chrN.maternal/paternal to SAMPLE#H#chrN" >> "$manifest.tmp"
    mv "$manifest.tmp" "$manifest"
done
