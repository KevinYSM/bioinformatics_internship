#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

params.fastq_directory="/data/local/proj/bioinformatics_project/data/raw/exome/*.fastq.gz"

params.outdir="/data/local/proj/bioinformatics_project/data/processed/qc_workflow_nextflow"

process FASTQC{
    input:
        file FASTQ_file
    output:
        path '*'

    """
    fastqc ${FASTQ_file}
    """
}

process MULTI{

    publishDir "${params.outdir}", mode: 'symlink'
    input: 
        file("*") 
    
    output:
        "*"
    
    """
    multiqc .
    """

    }

workflow{
    fastq_ch=Channel.fromPath(params.fastq_directory)
    multiqc_channel=FASTQC(fastq_ch).collect()
    multiqc_channel.view()
    MULTI(multiqc_channel)
    

}