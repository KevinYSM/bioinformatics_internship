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
samtools sort -n -m 2G -O BAM -o $(dirname "$bam")/$(basename "$bam" ".bam").sorted_by_name.bam "$bam"
