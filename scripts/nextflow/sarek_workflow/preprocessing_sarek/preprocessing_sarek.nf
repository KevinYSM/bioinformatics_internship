nextflow.enable.dsl=2
params.EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/unfinished_preprocessing_umis/*_R{1,2}_*.fastq.gz"
params.TEST_EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/50*_R{1,2}_*.fastq.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/"
params.disambiguate_outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/disambiguate/"
params.trimmed_dir="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis"
params.trimmed_umis="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"


process BWA_align_human {
    publishDir "${params.outdir}/bam", mode: 'symlink'
    maxForks 2
    input:
        file EXOME_trimmed_read_pair
    output:
        path "*.bam"
    """
    filename_human=\$(basename ${EXOME_trimmed_read_pair[1]} _Cut_0.fastq.gz)
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_human} ${EXOME_trimmed_read_pair[1]} ${EXOME_trimmed_read_pair[2]} | samtools sort -@7 -n -O BAM -o \${filename_human::-22}_human.sorted_by_name.bam -

    """
}


process BWA_align_mouse {
    publishDir "${params.outdir}/bam", mode: 'symlink'
    maxForks 2
    input:
        file EXOME_trimmed_read_pair
    output:
        file "*.bam"

    //https://biology.stackexchange.com/questions/59493/how-to-convert-bwa-mem-output-to-bam-format-without-saving-sam-file
    """
    filename_mouse=\$(basename ${EXOME_trimmed_read_pair[1]} _Cut_0.fastq.gz)
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_mouse} ${EXOME_trimmed_read_pair[1]} ${EXOME_trimmed_read_pair[2]} | samtools sort -@7 -n -O BAM -o \${filename_mouse::-22}_mouse.sorted_by_name.bam -
    
    """
}


process SAM_sort_name {
    maxForks 9
    publishDir "${params.outdir}/sam_sorted", mode: 'symlink'
    input:
        file ALIGNED_bam_file
    output:
        file '*.bam'
    """
    SAM_sort_name.sh $ALIGNED_bam_file
    """
}



process NGS_disambiguate {
        maxForks 20
    publishDir "${params.outdir}/disambiguate" , mode: 'copy'
    input:
        tuple val(key), path(BAM_files)
    output:
        file '*.disambiguatedSpeciesA.bam'
        file '*.txt'
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

//NEED TO SORT BY READ NAME HERE

process SAM_sort{
        maxForks 2
        publishDir "${params.outdir}/sam_sort_name_2" , mode: 'symlink'
        input:
                file DISAMBIGUATE_bam_file
                file TXT_file

        output:
                file '*.bam'

        """
        SAM_sort_name.sh $DISAMBIGUATE_bam_file
        """
}
process SAM_bam_to_fastq {
        maxForks 2
        publishDir "${params.outdir}/fastq" , mode: 'copy'
        input:
                file SORTED_bam_file
        output:
                path "*.fastq.gz"
        """
        samtools fastq $SORTED_bam_file -1 \$(basename "$SORTED_bam_file" .bam).read1.fastq.gz -2 \$(basename "$SORTED_bam_file" .bam).read2.fastq.gz  -n

        """
}





workflow{
    trimmed_read_pairs_ch=Channel.fromFilePairs(params.EXOME_trimmed_reads_directory, flat:true)
  
    DISAMBIGUATE_human=BWA_align_human(trimmed_read_pairs_ch)
    DISAMBIGUATE_mouse=BWA_align_mouse(trimmed_read_pairs_ch)
    DISAMBIGUATE_all=DISAMBIGUATE_human.mix(DISAMBIGUATE_mouse)
    DISAMBIGUATE_all.view()
    //DISAMBIGUATE_ch.view()

    
   
   
    
    DISAMBIGUATE_ch=DISAMBIGUATE_all.map{it ->[it.name.split('_')[0],it] }.groupTuple()
    DISAMBIGUATE_ch.view()
    SAM_sort_disambiguate_ch=NGS_disambiguate(DISAMBIGUATE_ch)
    BAM_to_fastq_ch=SAM_sort(SAM_sort_disambiguate_ch)
    SAM_bam_to_fastq(BAM_to_fastq_ch)
 
    //DISAMBIGUATE_human=DISAMBIGUATE_all.filter{file -> file =~/human.sam/}
    //DISAMBIGUATE_mouse=DISAMBIGUATE_all.filter{file -> file =~/mouse.sam/}
    //DISAMBIGUATE_human.view()
    //DISAMBIGUATE_mouse.view()
    //DISAMBIGUATE_ch=DISAMBIGUATE_all.map{it ->[it.name.split('_')[0],it] }.groupTuple()
    //DISAMBIGUATE_ch.view()
    //DISAMBIGUATE_mouse=SAM_sort_name_2(Channel.fromPath("/data/local/proj/bioinformatics_project/data/interim/exome/sam/batch1/*mouse.sam"))
  
    
    //DISAMBIGUATE_human=SAM_sort_name(STAR_ch_human)
    //DISAMBIGUATE_mouse=SAM_sort_name_2(STAR_ch_mouse)
    //DISAMBIGUATE_ch=DISAMBIGUATE_human.join(DISAMBIGUATE_mouse)
    //DISAMBIGUATE_ch=Channel.fromFilePairs("/data/local/proj/bioinformatics_project/data/interim/exome/sam/batch1/*{human,_mouse}.sam", flat: true)
    //DISAMBIGUATE_ch.view()


    //convert_to_fastq_ch=NGS_disambiguate(DISAMBIGUATE_ch)
    //SAM_bam_to_fastq(convert_to_fastq_ch)
    //SAM_sort_ch=NGS_disambiguate(DISAMBIGUATE_ch)
}

