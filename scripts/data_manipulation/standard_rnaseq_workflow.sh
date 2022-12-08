#!/usr/bin/env bash


#This scripts runs the standard RNA-seq preprocessing workflow on all fastq files within a directory
# STAR -> samtools -> HTSeq

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

PROJECT_DIRECTORY=/data/local/proj/bioinformatics_project

#STEP ZERO: Ensure file system is set-up correctly, generate log file
bash "$PROJECT_DIRECTORY"/scripts/tools/initdir_standard_rnaseq.sh

if [ ! -d "logs" ]
then
    mkdir logs
fi
1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}').log 2>&1 #how does this line of code work?

#STEP ONE: Perform STAR alignment on all read pairs
in_dir=$1
out_dir=/data/local/data/interim/rna/aligned

for R1 in $(find "$1"/ -name "*_R1_001.fastq.gz")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/standard_rnaseq_scripts/STAR_align_human.sh $R1
done

wait


#STEP TWO: Use Samtools to sort aligned .bam files by name
out_dir_2=/data/local/data/interim/rna/sorted_name
for ALIGNED in $(find "$out_dir"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/standard_rnaseq_scripts/SAM_sort_name.sh $ALIGNED $out_dir_2
done

wait


#STEP THREE: Use HtSeq to count read frequencies for genes.
out_dir_3=/data/loca/data/processed/rna/counts
for SORTED in $(find "$out_dir_2"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/standard_rnaseq_scripts/HTSEQ_count.sh $ALIGNED $out_dir_3
done
