#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/variant_calling_workflow_nextflow/"
//params.recal_bam_files="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/PCB*"
params.recal_bam_files="/data/local/proj/bioinformatics_project/data/interim/rna/base_recalibrated/*.recal.pass2.bam"

//params.RNA_aligned_directory="/data/local/proj/bioinformatics_project/data/interim/rna/aligned"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"


process GATK_haplotype_caller{
       
        input:
                file RECAL_bam_file
        output:
        path '*.vcf.gz'
        path '*.tbi'

        """
        variant_calling.sh ${RECAL_bam_file} ${params.fasta_human} 
        """
}

process GATK_variant_filtration{
        publishDir "${params.outdir}", mode: 'symlink'
        input:
                file VCF_file
                file TBI_file
        output:
        path '*.vcf.gz'

        """
        filter_variants.sh ${VCF_file} ${params.fasta_human} 
        """
}


workflow{
        GATK_haplotype_caller_ch=Channel.fromPath(params.recal_bam_files)
        GATK_variant_filtration_ch=GATK_haplotype_caller(GATK_haplotype_caller_ch)
        GATK_variant_filtration(GATK_variant_filtration_ch)

}

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}

