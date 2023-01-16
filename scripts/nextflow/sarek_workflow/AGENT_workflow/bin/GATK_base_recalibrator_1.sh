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
    -O $(basename "$bam" ".bam").recal.pass1.table
