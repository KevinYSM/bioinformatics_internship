#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

recal_pass2_table=$1

bam=$(basename "$recal_pass2_table" ".recal.pass2.table")".bam"

tmpdir="/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/""$(basename "$bam" ".bam")"
fasta=$2
outdir="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/output"

k1=$3
k2=$4


gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I "$outdir"/$(basename "$bam" ".bam").recal.pass1.bam \
    -R "$fasta" \
    --bqsr-recal-file $recal_pass2_table \
    -O $(basename "$bam" ".bam").recal.pass2.bam

if [ -f "$outdir"/$(basename "$bam" ".bam").recal.pass2.bam ]
then
    rm "$outdir"/$(basename "$bam" ".bam").recal.pass1.bam
    rm "$outdir"/$(basename "$bam" ".bam").recal.pass1.bam.bai
fi


rm -r "$tmpdir"