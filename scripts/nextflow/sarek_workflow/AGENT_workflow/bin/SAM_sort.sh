#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o verbose
bam=$1
samtools sort -O BAM -o $(basename "$bam" '.bam').sorted.bam "$bam" 
samtools index $(basename "$bam" ".bam").sorted.bam