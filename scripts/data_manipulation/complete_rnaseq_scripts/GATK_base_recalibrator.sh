#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1
output_dir=$2

tmpdir=tmp_split_$output_dir
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"

fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta

k1=/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz \
k2=/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I "$bam" \
    -R "$fasta" \
    --known-sites "$k1" \
    --known-sites "$k2" \
    -O $output_dir/$(basename "$bam" ".bam").recal_pass1.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I "$bam" \
    -R "$fasta" \
    --bqsr-recal-file "$bam".recal_pass1.table \
    -O $output_dir/$(basename "$bam" ".bam").recal.pass1.bam

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R "$fasta" \
    --known-sites "$k1" \
    --known-sites "$k2" \
    -O $output_dir/$(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R "$fasta" \
    --bqsr-recal-file $output_dir/$(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table \
    -O $output_dir/$(basename "$bam" ".bam").recal.pass2.bam

if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass2.bam ]
then
    rm $output_dir/$(basename "$bam" ".bam").recal.pass1.bam
    rm $output_dir/$(basename "$bam" ".bam").recal.pass1.bam.bai
fi


rm -r "$tmpdir"