#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input.bam> <k> <output.bam>" >&2
    exit 1
fi

INPUT_BAM="$1"
K="$2"
OUTPUT_BAM="$3"

if ! [[ "$K" =~ ^[0-9]+$ ]] || [[ "$K" -le 0 ]]; then
    echo "ERROR: k must be a positive integer." >&2
    exit 1
fi

samtools view -h "$INPUT_BAM" \
  | awk -v k="$K" '
      BEGIN { removed=0 }

      /^@/ {
          print
          next
      }

      {
          flag=$2

          if ((and(flag,4)==0) && removed<k) {
              removed++
              next
          }

          print
      }

      END {
          print "removed_mapped_records=" removed > "/dev/stderr"

          if (removed != k) {
              print "ERROR: requested " k \
                    " mapped records but removed " removed > "/dev/stderr"
              exit 2
          }
      }
    ' \
  | samtools view -b -o "$OUTPUT_BAM" -

samtools quickcheck -v "$OUTPUT_BAM"
