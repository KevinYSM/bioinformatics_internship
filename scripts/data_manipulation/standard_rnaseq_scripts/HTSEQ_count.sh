#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam="$1"
output_dir="$2"
#gtf=genome/Homo_sapiens.GRCh37.75.sorted.gtf
gtf=/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf

htseq-count -r name -q -f bam -s reverse -m intersection-strict -o $output_dir/$(basename "$bam" ".bam").s_reverse.gene_counts  
wait
