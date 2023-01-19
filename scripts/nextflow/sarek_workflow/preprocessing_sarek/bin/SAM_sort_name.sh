#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o verbose

bam=$1

samtools sort -n -@7 -O BAM -o $(basename "$bam" ".bam").sorted_by_name.bam "$bam"