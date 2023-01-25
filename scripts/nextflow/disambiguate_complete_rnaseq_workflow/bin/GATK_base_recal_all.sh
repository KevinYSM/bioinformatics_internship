#!/usr/bin/env bash
bam=$1
fasta=$2
dbsnp=$3
known_indels=$4


tmpdir=/data/local/proj/bioinformatics_project/scripts/nextflow/disambiguate_complete_rnaseq_workflow/temp_dir/$(basename "$bam" ".bam")
if [ -d "$tmpdir" ]
then
    rm -r "$tmpdir"
fi
mkdir -p "$tmpdir"


gatk BaseRecalibrator --tmp-dir "$tmpdir"\
        -I "$bam" \
        -R "$fasta" \
        --known-sites "$dbsnp" \
        --known-sites "$known_indels" \
        -O $(basename "$bam" ".bam").recal.pass1.table

 

    gatk ApplyBQSR --tmp-dir "$tmpdir"\
        -I "$bam" \
        -R "$fasta" \
        --bqsr-recal-file $(basename "$bam" ".bam").recal.pass1.table \
        -O $(basename "$bam" ".bam").recal.pass1.bam

 

    gatk BaseRecalibrator --tmp-dir "$tmpdir"\
        -I $(basename "$bam" ".bam").recal.pass1.bam \
        -R "$fasta" \
        --known-sites "$dbsnp" \
        --known-sites "$known_indels" \
        -O $(basename "$bam" ".bam").recal.pass2.table

 

    gatk ApplyBQSR --tmp-dir "$tmpdir"\
        -I $(basename "$bam" ".bam").recal.pass1.bam \
        -R "$fasta" \
        --bqsr-recal-file $(basename "$bam" ".bam").recal.pass2.table \
        -O $(basename "$bam" ".bam").recal.pass2.bam


rm $(basename "$bam" ".bam").recal.pass1.table
rm $(basename "$bam" ".bam").recal.pass1.bam

