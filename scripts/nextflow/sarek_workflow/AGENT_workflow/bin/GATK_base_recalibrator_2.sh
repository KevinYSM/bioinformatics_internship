#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

recal_pass1_bam=$1
bam=$(basename "$recal_pass1_bam" ".recal.pass1.bam")".bam"
tmpdir="/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/""$(basename "$bam" ".bam")"


fasta=$2

k1=$3
k2=$4

gatk BaseRecalibrator --tmp-dir "$tmpdir" \
    -I $recal_pass1_bam  \
    -R "$fasta" \
    --known-sites "$k1" \
    --known-sites "$k2" \
    -O $(basename "$bam" ".bam").recal.pass2.table