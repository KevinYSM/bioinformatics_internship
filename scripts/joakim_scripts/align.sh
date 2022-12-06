#!/bin/bash

input_fastq=$1
output_bam=$(basename $input_fastq ".fastq").bam

STAR -p 12 -t 2 $input_fastq -o ./results/$output_bam
