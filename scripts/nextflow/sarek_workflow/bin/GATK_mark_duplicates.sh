#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace
set -o verbose

bam=$1


gatk MarkDuplicates \
    --INPUT "$bam" \
    --OUTPUT $(basename "$bam" ".bam").MarkDuplicates.bam \
    --METRICS_FILE $(basename "$bam" ".bam").MarkDuplicates.metrics.txt
