#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1
output_dir=$2

samtools sort -O BAM -o $output_dir/$(basename "$bam" ".bam").sorted.bam "$bam"
samtools index $output_dir/$(basename "$bam" ".bam").sorted.bam(dirname "$bam")