params.EXOME_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/exome/*_R{1,2}_*.fastq.gz"
params.TEST_EXOME_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/exome/16-PDX*_R{1,2}_*.fastq.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/"
params.trimmed_dir="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"

process AGENT_trim_umis{ 
    maxForks 4
    
    input:
        file EXOME_read_pair
    output:
      
    """
    agent.sh ${EXOME_read_pair[1]} ${EXOME_read_pair[2]} ${params.trimmed_dir}
    """
}

process BWA_align_human{
    input:
        file EXOME_read_pair
    output:
        path "*.sam"
    """
    bwa mem -t 3 -K 100000000 -Y ${params.fast_human} read1.fastq.gz read2.fastq.gz
    """
}


process BWA_align_mouse{
    input:
        file EXOME_read_pair
    output:
        path "*.sam"
    """
    bwa mem -t 3 -K 100000000 -Y ${params.fasta_mouse} read1.fastq.gz read2.fastq.gz
    """
}

process SAM_sort_name{
    input:
        file ALIGNED_sam_file
    output:
        path '*.bam'
    """
    SAM_sort_name.sh $ALIGNED_sam_file
    """
}

process SAM_sort_name_2{
    input:
        file ALIGNED_sam_file
    output:
        path '*.bam'
    """
    SAM_sort_name.sh $ALIGNED_sam_file
    """
}

process NGS_disambiguate{
    input:
        path BAM_files
    output:
        path '*.disambiguatedSpeciesA.bam'
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

process convert_bam_to_fastq{
    input:
        file DISAMBIGUATED_bam_file
    output:
        path "*.fastq"
    """
    samtools fastq $DISAMBIGUATED_bam_file > \$(basename "$DISAMBIGUATED_bam_file" ".bam").fastq
    """
}

process SAREK_pipeline{
    input:
        file
    output:
        path
    """
    """
}

workflow{
    read_pairs_ch = Channel.fromFilePairs(params.TEST_EXOME_raw_reads_directory, flat: true)
    read_pairs_ch.view()
    BWA_align_ch=AGENT_trim_umis(read_pairs_ch)
    BWA_align_ch.view()

    //BWA_align_human(BWA_align_ch)
    //BWA_align_mouse(BWA_align_ch)
}