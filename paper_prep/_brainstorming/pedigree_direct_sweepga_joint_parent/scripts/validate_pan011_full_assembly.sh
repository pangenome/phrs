#!/usr/bin/env bash
set -euo pipefail

FASTA=${1:-/moosefs/erikg/phrs/recovery/prereq-restore-readable/PAN011.fa.gz}

gzip -t "$FASTA"
test -s "$FASTA.fai"
head "$FASTA.fai"

grep -q '^PAN011#1#chr1[[:space:]]' "$FASTA.fai"
grep -q '^PAN011#2#chr1[[:space:]]' "$FASTA.fai"
