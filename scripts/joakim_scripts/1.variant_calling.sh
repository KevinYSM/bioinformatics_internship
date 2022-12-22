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
fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta
#fasta=/home/joakim/store1/pinbot/data1/reference/GATK_resource_bundle/b37/human_g1k_v37.fasta

#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
#https://gatkforums.broadinstitute.org/gatk/discussion/8692/version-highlights-for-gatk-version-3-7
gatk HaplotypeCaller -R="$fasta" -I="$bam" --dont-use-soft-clipped-bases=true --standard-min-confidence-threshold-for-calling=20.0 -O=$(dirname "$bam")/$(basename "$bam").haplotypecaller.vcf.gz
