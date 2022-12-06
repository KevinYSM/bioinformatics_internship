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
#genome="$2"
#genome=/home/joakim/store1/pinbot/data1/reference/GATK_resource_bundle/b37/human_g1k_v37.fasta
genome=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta

tmpdir=tmp_split_$(basename "$bam")
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"

#https://gatkforums.broadinstitute.org/gatk/discussion/10800/gatk4-how-to-reassign-star-mapping-quality-from-255-to-60-with-splitncigarreads
#https://github.com/bcbio/bcbio-nextgen/issues/2080
gatk SplitNCigarReads --tmp-dir "$tmpdir" -R="$genome" -I="$bam" -O=$(dirname "$bam")/$(basename "$bam" ".bam").split.bam
#if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").split.bam ]
#then
#    rm "$bam"
#    rm "$bam".bai
#fi
rm -r "$tmpdir"
