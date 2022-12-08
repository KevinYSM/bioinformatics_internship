#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1
output_dir=$2

gatk MarkDuplicates \
    --INPUT "$bam" \
    --OUTPUT $output_dir/$(basename "$bam" ".bam").MarkDuplicates.bam \
    --METRICS_FILE $output_dir/$(basename "$bam" ".bam").MarkDuplicates.metrics.txt
samtools index -@ 5 $output_dir/$(basename "$bam" ".bam").MarkDuplicates.bam