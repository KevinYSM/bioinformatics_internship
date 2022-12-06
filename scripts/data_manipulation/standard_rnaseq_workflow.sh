#!/usr/bin/env bash


#This scripts runs the standard RNA-seq preprocessing workflow on all fastq files within a directory
# STAR -> samtools -> HTSeq

set -u
set -o errexit
#set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

#STEP ZERO: Ensure file system is set-up correctly


#STEP ONE: Perform STAR alignment on all read pairs

#STEP TWO: Use Samtools to sort aligned .bam files by name

#STEP THREE: Use HtSeq to count read frequencies for genes
