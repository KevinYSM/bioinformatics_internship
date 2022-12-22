#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1

tmpdir=/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/$(basename "$bam" ".bam")
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"

fasta=$2

k1=$3
k2=$4

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I "$bam" \
    -R "$fasta" \
    --known-sites "$k1" \
    --known-sites "$k2" \
    -O $(basename "$bam" ".bam").recal_pass1.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I "$bam" \
    -R "$fasta" \
    --bqsr-recal-file "$bam".recal_pass1.table \
    -O $(basename "$bam" ".bam").recal.pass1.bam

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R "$fasta" \
    --known-sites "$k1" \
    --known-sites "$k2" \
    -O $(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table

#--known-sites="$k3" \

gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass1.bam \
    -R "$fasta" \
    --bqsr-recal-file $output_dir/$(basename "$bam" ".bam").recal.pass1.bam.recal_pass2.table \
    -O $(basename "$bam" ".bam").recal.pass2.bam

if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").recal.pass2.bam ]
then
    rm $output_dir/$(basename "$bam" ".bam").recal.pass1.bam
    rm $output_dir/$(basename "$bam" ".bam").recal.pass1.bam.bai
fi


rm -r "$tmpdir"