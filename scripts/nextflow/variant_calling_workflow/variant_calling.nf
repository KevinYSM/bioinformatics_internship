#!/usr/bin/env nextflow
// Declare syntax version

//Used to produce .vcf files following disambiguate workflow
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/variant_calling_workflow_nextflow/"
//params.recal_bam_files="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/PCB*"
params.recal_bam_files="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/*.recal.pass2.bam"

//params.RNA_aligned_directory="/data/local/proj/bioinformatics_project/data/interim/rna/aligned"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"


process GATK_haplotype_caller{
        publishDir "${params.outdir}/GATK_haplotype_caller", mode: 'symlink'
        input:
                path RECAL_bam_file
        output:
                path '*.vcf.gz'
                path '*.tbi'

        """
        variant_calling.sh ${RECAL_bam_file} ${params.fasta_human} 
        """
}

process GATK_variant_filtration{
        publishDir "${params.outdir}/GATK_variant_filtration", mode: 'symlink'
        input:
                file VCF_file
                file TBI_file
        output:
                path '*.vcf.gz'

        """
        filter_variants.sh ${VCF_file} ${params.fasta_human} 
        """
}

process VEP {
    publishDir "${params.outdir}/VEP", mode: 'symlink'

    input:
        path vcf_gz
    output:
        path '*.vcf'

    script:
    """
    vcf=\$(basename "${vcf_gz}" ".vcf.gz").vcf
    bgzip -d -c -f "${vcf_gz}" > "\$vcf"

 

    if [ ! -d output_vep_updated ]
    then
        mkdir output_vep_updated
    fi
    singularity exec \
        -B \$(pwd)/output_vep_updated:/output_vep_updated \
        -B /data/vep_cache:/.vep \
        -B "\$vcf":/\$(basename "\$vcf") \
        /home/ubuntu/vep.sif /opt/vep/src/ensembl-vep/vep \
        --species homo_sapiens \
        --assembly GRCh38 \
        --offline \
        --cache \
        --dir /.vep \
        --input_file /\$(basename "\$vcf") \
        --output_file /output_vep_updated/\$(basename "\$vcf" ".vcf").ann.vcf \
        --everything \
        --vcf \
        --fasta /.vep/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --force_overwrite
    """
}

 

process VCF2MAF {
    publishDir "${params.outdir}/MAF", mode: 'copy'

    input:
        path vcf
    output:
        path("*.maf")
    script:
    """
    tumor_id=\$(basename "${vcf}" ".ann.vcf" |  awk -F"_L004" '{print \$1}')
    vcf2maf.pl \
        --inhibit-vep \
        --input-vcf "${vcf}" \
        --output-maf \$(basename "${vcf}" ".vcf").maf \
        --tumor-id "\$tumor_id" \
        --ref-fasta /data/vep_cache/homo_sapiens/105_GRCh38/Homo_sapiens_assembly38.fasta \
        --ncbi-build GRCh38
    """
}

workflow{
        GATK_haplotype_caller_ch=Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/*.recal.pass2.bam")
        GATK_variant_filtration_ch=GATK_haplotype_caller(GATK_haplotype_caller_ch)
        VEP_ch=GATK_variant_filtration(GATK_variant_filtration_ch)
        VCF_ch=VEP(VEP_ch)
       
        VCF2MAF(VCF_ch)

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

