#!/usr/bin/env bash

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

ulimit -n 65535 #temporary file limit

n_cores=30
fasta=/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
gtf=/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
genomeDir=/data/local/proj/bioinformatics_project/data/external/genome
out_dir=/data/local/data/interim/rna/aligned


sample=$(basename "$1" | awk -F"_R1" '{print $1}')
R2=$(echo "$1" | sed s/_R1_/_R2_/g)
ID=$(zcat "$1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $3}').$(zcat "$1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $4}')
#ID=$(zcat "$R1" | head -n 1 | awk -F":" '{print $1}' | awk -F"@" '{print $2}').$(zcat "$R1" | head -n 1 | awk -F":" '{print $2}')
SM="$sample"
PL=ILLUMINA
LB="$SM" # Assume that each sample has only been subjected to one library preparation, which has been used for both runs (no other information has been given to indicate otherwise).

maxReadLength=$(awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(gunzip -c "$1") | awk '{print $1}' | tail -n 1)

STAR --runThreadN "$n_cores" --outBAMsortingThreadN "$n_cores" --twopassMode Basic --genomeDir "$genomeDir" --outFilterType BySJout --readFilesIn "$1" "$R2" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$out_dir"/"$sample"."$ID". --outSAMattrRGline ID:${ID} PL:${PL} LB:${LB} SM:${SM} --sjdbOverhang $(expr $maxReadLength - 1) --sjdbGTFfile "$gtf" --outSAMmapqUnique 60 --readFilesCommand zcat

