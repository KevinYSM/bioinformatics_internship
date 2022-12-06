#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}')."$RANDOM".log 2>&1

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam="$1"

#tmpdir=tmp_mkdup_$(basename "$bam")
#if [ -d "$tmpdir" ]
#then
#    rm -r "$tmpdir"
#fi
#mkdir -p "$tmpdir"

gatk MarkDuplicates \
    --INPUT="$bam" \
    --OUTPUT=$(dirname "$bam")/$(basename "$bam" ".bam").MarkDuplicates.bam \
    --METRICS_FILE=$(dirname "$bam")/$(basename "$bam" ".bam").MarkDuplicates.metrics.txt
samtools index -@ 5 $(dirname "$bam")/$(basename "$bam" ".bam").MarkDuplicates.bam

#rm -r "$tmpdir"
