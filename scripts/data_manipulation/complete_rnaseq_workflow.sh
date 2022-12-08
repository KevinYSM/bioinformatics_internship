#!/usr/bin/env bash

#This scripts runs the complete RNA-seq preprocessing workflow on all fastq files within a directory
# STAR -> samtools (index)-> gatk MarkDuplicates -> gatk SplitNCigarReads -> gatk Recalibration

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

PROJECT_DIRECTORY=/data/local/proj/bioinformatics_project

#STEP ZERO: Ensure file system is set-up correctly, generate log file
bash "$PROJECT_DIRECTORY"/scripts/tools/initdir_complete_rnaseq.sh

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

#STEP TWO: Sort .bam files by index using samtools

out_dir_2=/data/local/data/interim/rna/sorted_index
for ALIGNED in $(find "$out_dir"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/complete_rnaseq_scripts/SAM_sort_index.sh $ALIGNED $out_dir_2
done

wait

#STEP THREE:  Mark duplicates using gatk

out_dir_3=/data/local/data/interim/rna/mark_duplicates
for SORTED in $(find "$out_dir_2"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/complete_rnaseq_scripts/GATK_mark_duplicates.sh $SORTED $out_dir_3
done

wait

#STEP FOUR: gatk SplitNCigarReads??

out_dir_4=/data/local/data/interim/rna/splitNtrim

for MARK_DUPLICATES in $(find "$out_dir_3"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/complete_rnaseq_scripts/GATK_split_N_cigar_reads.sh
done

wait
#STEP FIVE: gatk recalibration??

out_dir_5=/data/local/data/processed/rna/recalibrated

for SPLIT_N_TRIM in $(find "$out_dir_4"/ -name "*.bam")
do
        bash "$PROJECT_DIRECTORY"/scripts/data_manipulation/complete_rnaseq_scripts/GATK_base_recalibration.sh
done