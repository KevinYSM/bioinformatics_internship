params.EXOME_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/exome/*_R{1,2}_*.fastq.gz"
params.TEST_EXOME_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/exome/16-PDX*_R{1,2}_*.fastq.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/"
params.trimmed_dir="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis"
params.trimmed_umis="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"
//bwa mem -t 12 -K 100000000 -Y /data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna  16-PDX_S1_L001_R1_001.1673495818136_Cut_0.fastq.gz 16-PDX_S1_L001_R2_001.1673495818136_Cut_0.fastq.gz > /data/local/proj/bioinformatics_project/data/interim/exome/bwa_mem/16-PDX_S1_L001_cut_mouse.sam
process AGENT_trim_umis{ 
    maxForks 5
    
    input:
        file EXOME_read_pair
    output:
      
    """
    agent.sh ${EXOME_read_pair[1]} ${EXOME_read_pair[2]} ${params.trimmed_dir}
    """
}



workflow{
    //read_pairs_ch = Channel.fromFilePairs(params.EXOME_raw_reads_directory, flat: true)
    //BWA_align_ch=AGENT_trim_umis(read_pairs_ch)
    trimmed_read_pairs_ch=Channel.fromFilePairs(params.trimmed_umis, flat:true)
    trimmed_read_pairs_ch.view()
    
}