#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}').log 2>&1

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam="$1"
samtools sort -O BAM -o $(dirname "$bam")/$(basename "$bam" ".bam").sorted.bam "$bam"
samtools index $(dirname "$bam")/$(basename "$bam" ".bam").sorted.bam
#if [ -f "$(dirname "$bam")/$(basename "$bam" ".bam").sorted.bam.bai" ]
#then
#    rm "$bam"
#fi
