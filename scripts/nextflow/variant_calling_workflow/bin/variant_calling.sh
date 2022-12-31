#!/usr/bin/env bash


bam="$1"
fasta=$2
#fasta=/home/joakim/store1/pinbot/data1/reference/GATK_resource_bundle/b37/human_g1k_v37.fasta

#https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
#https://gatkforums.broadinstitute.org/gatk/discussion/8692/version-highlights-for-gatk-version-3-7
gatk HaplotypeCaller -R "$fasta" -I "$bam" --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20.0 -O $(basename "$bam").haplotypecaller.vcf.gz
