#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam=$1
output_dir=$2

fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta

tmpdir=tmp_split_$output_dir
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"

#https://gatkforums.broadinstitute.org/gatk/discussion/10800/gatk4-how-to-reassign-star-mapping-quality-from-255-to-60-with-splitncigarreads
#https://github.com/bcbio/bcbio-nextgen/issues/2080
gatk SplitNCigarReads --tmp-dir "$tmpdir" -R "$fasta" -I "$bam" -O $output_dir/$(basename "$bam" ".bam").split.bam
#if [ -f $(dirname "$bam")/$(basename "$bam" ".bam").split.bam ]
#then
#    rm "$bam"
#    rm "$bam".bai
#fi
rm -r "$tmpdir"
