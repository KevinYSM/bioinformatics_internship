#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose
ulimit -n 65535 #temporary file limit
bqsr=$1
outdir="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/output"

bam="$outdir"/$(basename "$bqsr" ".recal.pass1.table")".bam"
tmpdir=/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/$(basename "$bam" ".bam")
fasta=$2
k1=$3
k2=$4



gatk ApplyBQSR --tmp-dir "$tmpdir" \
    -I "$bam" \
    -R "$fasta" \
    --bqsr-recal-file "$bqsr" \
    -O $(basename "$bam" ".bam").recal.pass1.bam