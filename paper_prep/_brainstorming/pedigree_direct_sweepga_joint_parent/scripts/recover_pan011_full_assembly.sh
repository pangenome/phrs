#!/usr/bin/env bash
set -euo pipefail

OUT_DIR=${1:-/moosefs/erikg/phrs/recovery/prereq-restore-readable}
URL=${PAN011_PUBLIC_FASTA_URL:-https://public.gi.ucsc.edu/~mcechova/pedigree/assemblies/v1.1/assembly.v1.1.PAN011.diploid.fa}

mkdir -p "$OUT_DIR"

tmp="$OUT_DIR/PAN011.fa.gz.tmp"
out="$OUT_DIR/PAN011.fa.gz"

curl -L "$URL" |
  awk 'BEGIN{FS="[ .]"}
       /^>PAN011\.chr/ {
         chrom=$2
         hap=$3
         sub(/^haplotype/, "", hap)
         print ">PAN011#" hap "#" chrom
         next
       }
       {print}' |
  bgzip -@ "${BGZIP_THREADS:-4}" -c > "$tmp"

mv "$tmp" "$out"
samtools faidx "$out"

gzip -t "$out"
head "$out.fai"
