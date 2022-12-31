#!/usr/bin/env bash
vcf=$1
fasta=$2
#fasta="$2"

gatk VariantFiltration -R "$fasta" -V "$vcf" --cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" -O $(basename "$vcf" ".vcf.gz").filtered.vcf.gz
