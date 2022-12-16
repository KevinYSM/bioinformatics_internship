#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o verbose

bam=$1

samtools sort -n -O BAM -o $(basename "$bam" ".bam").sorted_by_name.bam "$bam"