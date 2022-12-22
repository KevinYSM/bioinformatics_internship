#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}')."$RANDOM".log 2>&1

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

vcf="$1"
fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta
#fasta="$2"

gatk VariantFiltration -R="$fasta" -V="$vcf" --cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" -O=$(dirname "$vcf")/$(basename "$vcf" ".vcf.gz").filtered.vcf.gz
