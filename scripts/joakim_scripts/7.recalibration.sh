#!/bin/bash
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

#fasta=~/store1/pinbot/data1/reference/GATK_resource_bundle/b37/human_g1k_v37.fasta
#k1=~/store1/pinbot/data1/reference/GATK_resource_bundle/b37/dbsnp_138.b37.vcf.gz
#k2=~/store1/pinbot/data1/reference/GATK_resource_bundle/b37/1000G_phase1.indels.b37.vcf.gz
#k3=~/store1/pinbot/data1/reference/GATK_resource_bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz

tmpdir=tmp_recal_$(basename "$bam")
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"


fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta

k1=/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz \
k2=/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I="$bam" \
    -R="$fasta" \
    --known-sites="$k1" \
    --known-sites="$k2" \
    -O="$bam".recal_pass1.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I="$bam" \
    -R="$fasta" \
    --bqsr-recal-file="$bam".recal_pass1.table \
    -O=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R="$fasta" \
    --known-sites="$k1" \
    --known-sites="$k2" \
    -O=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R="$fasta" \
    --bqsr-recal-file=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table \
    -O=$(dirname "$bam")/$(basename "$bam" ".bam").recal.pass2.bam

if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass2.bam ]
then
    rm $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam
    rm $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam.bai
fi
