#!/usr/bin/env bash

set -u
set -o errexit
set -o nounset
set -o xtrace
set -o verbose

ulimit -n 65535

n_cores=30

fasta=/data/local/reference/GATK_resource_bundle/hg38/hg38/Homo_sapiens_assembly38.fasta
gtf=/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
genomeDir=/data/local/proj/melanoma/Pipelines/rna/preprocessing/genome2

if [ ! -f "$genomeDir"/Genome ]
then #run once for a given genome
	STAR --runMode genomeGenerate --genomeDir "$genomeDir" --genomeFastaFiles "$fasta" --runThreadN "$n_cores"
fi
