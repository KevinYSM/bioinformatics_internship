params.vcf_files="/data/local/proj/bioinformatics_project/data/processed/variant_calling_workflow_nextflow/*.vcf.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/vcf_processing_workflow/"

process VEP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val(vcf_gz)
    output:
        path("output_vep_updated/*.ann.vcf")

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
    publishDir "${params.outdir}", mode: 'copy'

    input:
        val(vcf)
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
    VEP_ch=Channel.fromPath(params.vcf_files)
    VCF_ch=VEP(VEP_ch)
    VCF_ch.view()
    VCF2MAF(VCF_ch)
}