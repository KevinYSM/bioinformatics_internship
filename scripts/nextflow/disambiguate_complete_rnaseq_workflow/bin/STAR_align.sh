#!/usr/bin/env bash

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

ulimit -n 65535 #temporary file limit

n_cores=30

#out_dir=/data/local/data/interim/rna/aligned oops
#out_dir=/data/local/proj/bioinformatics_project/data/interim/rna/aligned


sample=$(basename "$1" | awk -F"_R1" '{print $1}')
R1=$1
R2=$2
fasta=$3
gtf=$4
genomeDir=$5
ID=$(zcat "$1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $3}').$(zcat "$1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $4}')
#ID=$(zcat "$R1" | head -n 1 | awk -F":" '{print $1}' | awk -F"@" '{print $2}').$(zcat "$R1" | head -n 1 | awk -F":" '{print $2}')
SM="$sample"
PL=ILLUMINA
LB="$SM" # Assume that each sample has only been subjected to one library preparation, which has been used for both runs (no other information has been given to indicate otherwise).

if [ ! -f "$genomeDir"/Genome ]
then #run once for a given genome
    STAR --runMode genomeGenerate --genomeDir "$genomeDir" --genomeFastaFiles "$fasta" --runThreadN "$n_cores"
fi

maxReadLength=$(awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(gunzip -c "$1") | awk '{print $1}' | tail -n 1)

STAR --runThreadN "$n_cores" --outBAMsortingThreadN "$n_cores" --twopassMode Basic --genomeDir "$genomeDir" --outFilterType BySJout --readFilesIn "$R1" "$R2" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$sample"."$ID". --outSAMattrRGline ID:${ID} PL:${PL} LB:${LB} SM:${SM} --sjdbOverhang $(expr $maxReadLength - 1) --sjdbGTFfile "$gtf" --outSAMmapqUnique 60 --readFilesCommand zcat

