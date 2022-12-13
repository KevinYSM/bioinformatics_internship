#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1

samtools index $(basename "$bam" ".bam").sorted.bam(dirname "$bam")