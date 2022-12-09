#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/standard_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/"

params.RNA_raw_reads_directory='/data/local/proj/bioinformatics_project/data/raw/rna/*_R{1,2}_*.fastq.gz'
params.fasta="data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.gtf="data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
params.genomeDir="data/local/proj/bioinformatics_project/data/external/genome"

process STAR_align_human{
        input:
        file RNA_read_pair
        

        output:
        path "*.bam"

        """
        STAR_align_human.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]}
        """

}

process SAM_sort_name{
        input:
        file ALIGNED_bam_file

        output:
        path '*.bam'

        """
        SAM_sort_name.sh $ALIGNED_bam_file
        """

}

process HTSEQ_count{
        publishDir params.outdir, mode: 'copy'
        input:
        file SORTED_bam_file

        output:
        path '*.gene_counts'

        """
        HTSEQ_count.sh $SORTED_bam_file
        """
}

process VIEW_1 {
        publishDir params.outdir, mode:'copy'
        input:
        file

        output:
        path '*.txt'

        """
        echo '$x'>'$x'.txt
        """
}

workflow{
        read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        STAR_ch=STAR_align_human(read_pairs_ch)
        SAM_ch=SAM_sort_name(STAR_ch)
        HTSEQ_count(SAM_ch)
        
}