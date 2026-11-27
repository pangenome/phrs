#!/bin/bash
#SBATCH -c 1
#SBATCH -t 00:10:00
# Can the frozen Pore-C rerun log (28,856,811,520 bytes since before 10:21 MST, job 2887442) still grow? Append one line
# from a compute node and report the exit code and error message of the write, and the size before and after.
# usage: sbatch append_test.sbatch FILE
set -uxo pipefail
F=$1; echo "before: $(stat -c %s $F) bytes"
echo "append test $(date -Is) from $(hostname)" >> $F; rc=$?
echo "write exit code $rc; after: $(stat -c %s $F) bytes"
