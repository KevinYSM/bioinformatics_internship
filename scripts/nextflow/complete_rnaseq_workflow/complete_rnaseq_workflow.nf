#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"


params.RNA_aligned_directory='data/local/proj/bioinformatics_project/data/interim/rna/aligned'
params.RNA_raw_reads_directory='/data/local/proj/bioinformatics_project/data/raw/rna/*_R{1,2}_*.fastq.gz'
params.fasta="data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.gtf="data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
params.genomeDir="data/local/proj/bioinformatics_project/data/external/genome"

params.k1="/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz"
params.k2="/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"


process SAM_sort{
        input:
        file ALIGNED_bam_file

        output:
        path '*.bam'

        """
      
        SAM_sort.sh $ALIGNED_bam_file
        """

}

process SAM_index{
        input:
        file SORTED_bam_file

        output:
        path '*.txt'

        """
        SAM_index.sh $SORTED_bam_file
        """
}

process GATK_mark_duplicates{
        input:
        file INDEXED_bam_file

        output:
        path '*.bam'

        """
        GATK_mark_duplicates.sh $INDEXED_bam_file 
        """
}

process GATK_split_ch{
        input:
        file DUPLICATES_bam_file

        output:
        path '*.bam'
        """
        GATK_split_N_cigar_reads.sh $DUPLICATES_bam_file ${params.fasta}
        """
}

process GATK_base_recalibration{
        publishDir params.outdir, mode:'copy'
        input:
        file SPLIT_bam_file

        output:
        path '*'

        """
        GATK_base_recalibrator.sh $SPLIT_bam_file ${params.fasta} ${k1} ${k2}
        """
}



workflow{
        def aligned_bam_files=Channel.fromPath('/data/local/proj/bioinformatics_project/data/interim/rna/aligned/*out.bam')
        SAM_index_ch=SAM_sort(aligned_bam_files)
        GATK_duplicates_ch=SAM_index(SAM_index_ch)
        GATK_split_ch=GATK_mark_duplicates(GATK_duplicates_ch)
        GATK_base_recalibration(GATK_split_ch)
        //read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        //STAR_ch=STAR_align_human(read_pairs_ch)
        //SAM_ch=SAM_sort_name(STAR_ch)
        //HTSEQ_count(SAM_ch)
        
}