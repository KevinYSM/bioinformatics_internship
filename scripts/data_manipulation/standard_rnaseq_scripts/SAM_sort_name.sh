#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1
output_dir=$2
samtools sort -n -m 10G -O BAM -o $output_dir/$(basename "$bam" ".bam").sorted_by_name.bam "$bam"