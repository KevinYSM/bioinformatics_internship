#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}').log 2>&1

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

ulimit -n 65535

n_cores=50
fasta=genome_mouse/GRCm38.primary_assembly.genome.fa
gtf=genome_mouse/gencode.vM22.annotation.gtf
genomeDir=$(dirname "$fasta")

if [ ! -f "$genomeDir"/Genome ]
then
    STAR --runMode genomeGenerate --genomeDir "$genomeDir" --genomeFastaFiles "$fasta" --runThreadN "$n_cores"
fi

in_dir=fastqdir/
out_dir=results/mouse/
if [ ! -d "$out_dir" ]
then
    mkdir -p "$out_dir"
fi
for R1 in $(find "$in_dir"/ -name "*_R1_001.fastq.gz")
do
    sample=$(basename "$R1" | awk -F"_R1" '{print $1}')
    R2=$(echo "$R1" | sed s/_R1_/_R2_/g)

    # Define read group
    ID=$(zcat "$R1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $3}').$(zcat "$R1" | head -n 1 | awk -F"@" '{print $2}' | awk -F":" '{print $4}')
    #ID=$(zcat "$R1" | head -n 1 | awk -F":" '{print $1}' | awk -F"@" '{print $2}').$(zcat "$R1" | head -n 1 | awk -F":" '{print $2}')
    SM="$sample"
    PL=ILLUMINA
    LB="$SM" # Assume that each sample has only been subjected to one library preparation, which has been used for both runs (no other information has been given to indicate otherwise).

    maxReadLength=$(awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' <(gunzip -c "$R1") | awk '{print $1}' | tail -n 1)

    STAR --runThreadN "$n_cores" --outBAMsortingThreadN "$n_cores" --twopassMode Basic --genomeDir "$genomeDir" --outFilterType BySJout --readFilesIn "$R1" "$R2" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$out_dir"/"$sample"."$ID". --outSAMattrRGline ID:${ID} PL:${PL} LB:${LB} SM:${SM} --sjdbOverhang $(expr $maxReadLength - 1) --sjdbGTFfile "$gtf" --outSAMmapqUnique 60 --readFilesCommand zcat
done
