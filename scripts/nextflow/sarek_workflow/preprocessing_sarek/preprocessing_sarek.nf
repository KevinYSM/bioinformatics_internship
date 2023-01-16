params.EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.TEST_EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/50*_R{1,2}_*.fastq.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/"
params.disambiguate_outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/disambiguate/"
params.trimmed_dir="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis"
params.trimmed_umis="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"
process BWA_align_human {
    maxForks 2
    input:
        file EXOME_trimmed_read_pair
    output:
      
    """
    filename_human=\$(basename ${EXOME_trimmed_read_pair[1]} _Cut_0.fastq.gz)
   
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_human} ${EXOME_trimmed_read_pair[1]} ${EXOME_trimmed_read_pair[2]} > /data/local/proj/bioinformatics_project/data/interim/exome/sam/\${filename_human::-21}_human.sam
    """
}


process BWA_align_mouse {
    maxForks 2
    input:
        file EXOME_trimmed_read_pair
    output:
        
    """
    filename_mouse=\$(basename ${EXOME_trimmed_read_pair[1]} _Cut_0.fastq.gz)}
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_mouse} ${EXOME_trimmed_read_pair[1]} ${EXOME_trimmed_read_pair[2]} > /data/local/proj/bioinformatics_project/data/interim/exome/sam/\${filename_mouse::-21}_mouse.sam
    """
}

process SAM_sort_name {
    input:
        file ALIGNED_sam_file
    output:
        path '*.bam'
    """
    SAM_sort_name.sh $ALIGNED_sam_file
    """
}

process SAM_sort_name_2 {
    input:
        file ALIGNED_sam_file
    output:
        path '*.bam'
    """
    SAM_sort_name.sh $ALIGNED_sam_file
    """
}

process NGS_disambiguate {
    publishDir "${params.disambiguate_outdir}" , mode: 'copy', pattern: "*"//save disambiguatedSpeciesA.bam
    input:
        path BAM_files
    output:
        path '*'
    """
       
        if [ ${BAM_files[1]} == *"human"* ]
        then
                NGS_disambiguate.sh ${BAM_files[1]} ${BAM_files[0]}
               
        elif [ ${BAM_files[1]} == *"mouse"* ]      
        then
                NGS_disambiguate.sh ${BAM_files[0]} ${BAM_files[1]}

        else
                echo "Error"
        fi

        """
}

process SAM_bam_to_fastq {
    publishDir "${params.outdir}" , mode: 'copy'
    input:
        file DISAMBIGUATED_bam_file
    output:
        path "*.fastq*"
    """
    samtools fastq $DISAMBIGUATED_bam_file > \$(basename "$DISAMBIGUATED_bam_file" ".bam").fastq
    """
}



process SAREK_pipeline {
    input:
        file
    output:
        path
    """
    """
}

workflow{
    trimmed_read_pairs_ch=Channel.fromFilePairs(params.EXOME_trimmed_reads_directory, flat:true)
    trimmed_read_pairs_ch.view();
    BWA_align_human(trimmed_read_pairs_ch)
    BWA_align_mouse(trimmed_read_pairs_ch)
    
  
    
    //DISAMBIGUATE_ch=DISAMBIGUATE_human.merge(DISAMBIGUATE_mouse)
    //DISAMBIGUATE_ch.view()
    //convert_to_fastq_ch=NGS_disambiguate(DISAMBIGUATE_ch).filter(~/.disambiguatedSpeciesA.bam/)
    //SAM_bam_to_fastq(convert_to_fastq_ch)
    //SAM_sort_ch=NGS_disambiguate(DISAMBIGUATE_ch)
}