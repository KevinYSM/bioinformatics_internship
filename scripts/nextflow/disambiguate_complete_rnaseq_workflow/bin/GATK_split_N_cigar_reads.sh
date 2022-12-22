#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1

fasta=$2

tmpdir=/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/$(basename "$bam" ".bam")
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"

#https://gatkforums.broadinstitute.org/gatk/discussion/10800/gatk4-how-to-reassign-star-mapping-quality-from-255-to-60-with-splitncigarreads
#https://github.com/bcbio/bcbio-nextgen/issues/2080
gatk SplitNCigarReads --tmp-dir "$tmpdir" -R "$fasta" -I "$bam" --output $(basename "$bam" ".bam").split.bam
#if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").split.bam ]
#then
#    rm "$bam"
#    rm "$bam".bai
#fi
rm -r "$tmpdir"
