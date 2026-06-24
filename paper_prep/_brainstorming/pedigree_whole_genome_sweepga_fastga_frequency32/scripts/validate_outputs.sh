#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

test -s "$PACKAGE_DIR/summaries/sweepga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/fastga_binary.tsv"
test -s "$PACKAGE_DIR/summaries/slurm_jobs.tsv"
test -s "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
test -s "$PACKAGE_DIR/summaries/input_manifest.tsv"

grep -R -- '--fastga-frequency 32' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '--scaffold-jump 0' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/dev/shm' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- '/home/erikg/.cargo/bin/sweepga' "$PACKAGE_DIR/logs" >/dev/null
grep -R -- 'FastGA -1:.* -f32 ' "$PACKAGE_DIR/logs" >/dev/null
grep -q $'PAN027_chr9q_chr3q_PHR_candidate\tPAN027pat_vs_PAN011_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"
grep -q $'PAN028_chr9q_chr3q_PHR_candidate\tPAN028mat_vs_PAN027_joint' "$PACKAGE_DIR/summaries/raw_chr3_support.tsv"

n_raw="$(find "$PACKAGE_DIR/raw_paf" -type f -name '*.paf.gz' | wc -l)"
if [[ "$n_raw" -ne 3 ]]; then
    echo "expected 3 raw PAFs, found $n_raw" >&2
    exit 1
fi

while IFS= read -r paf; do
    gzip -t "$paf"
    test -s "$paf.sha256"
    sha256sum -c "$paf.sha256"
    case "$(basename "$paf")" in
        *.frequency32_many_many_j0.paf.gz|*.sweepga_frequency32_many_many_j0.paf.gz) ;;
        *)
            echo "raw PAF name does not clearly include frequency32_many_many_j0: $paf" >&2
            exit 1
            ;;
    esac
done < <(find "$PACKAGE_DIR/raw_paf" -type f -name '*.paf.gz' | sort)

if git -C "$PACKAGE_DIR/../../.." status --short -- submission | grep -q .; then
    echo "submission/ was modified" >&2
    exit 1
fi

echo "Validated frequency32 sweepGA/FastGA raw whole-genome outputs"
