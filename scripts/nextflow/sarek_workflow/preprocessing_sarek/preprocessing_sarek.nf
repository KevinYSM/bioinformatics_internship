nextflow.enable.dsl=2
params.EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.TEST_EXOME_trimmed_reads_directory="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/50*_R{1,2}_*.fastq.gz"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/"
params.disambiguate_outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/disambiguate/"
params.trimmed_dir="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis"
params.trimmed_umis="/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/*_R{1,2}_*.fastq.gz"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"

params.SAM_completed_mouse=
[
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/36-PDX_S3_L001_R1_001.1673520115986_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/36-PDX_S3_L001_R2_001.1673520115986_Cut_0.fastq.gz" 
    ]
    ,
    [ "/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB136-TIL_S15_L002_R1_001.1673520035253_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB136-TIL_S15_L002_R2_001.1673520035253_Cut_0.fastq.gz"
    ]
    ,
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB36-TIL_S11_L002_R1_001.1673521649889_Cut_0.fastq.gz"
    , "/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB36-TIL_S11_L002_R2_001.1673521649889_Cut_0.fastq.gz"
    ]
]

params.SAM_completed_human=
[
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/39-PDX_S4_L001_R1_001.1673522893147_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/39-PDX_S4_L001_R2_001.1673522893147_Cut_0.fastq.gz"
    ]
    ,
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/54-PDX_S6_L001_R1_001.1673524805153_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/54-PDX_S6_L001_R2_001.1673524805153_Cut_0.fastq.gz"
    ]
    ,
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB-09-PDX_S1_L001_R1_001.1673518540275_Cut_0.fastq.gz",
    "/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB-09-PDX_S1_L001_R2_001.1673518540275_Cut_0.fastq.gz"
    ]
    ,
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB-30-PDX_S6_L001_R1_001.1673520055765_Cut_0.fastq.gz"
    , "/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB-30-PDX_S6_L001_R2_001.1673520055765_Cut_0.fastq.gz"
    ],
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB37-TIL_S13_L002_R1_001.1673522037838_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB37-TIL_S13_L002_R2_001.1673522037838_Cut_0.fastq.gz"
    ],
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB57-PDX_S3_L001_R1_001.1673523288426_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB57-PDX_S3_L001_R2_001.1673523288426_Cut_0.fastq.gz"
    ],
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB64-PDX_S1_L001_R1_001.1673524805265_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB64-PDX_S1_L001_R2_001.1673524805265_Cut_0.fastq.gz"
    ],
    ["/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB75-PDX_S5_L001_R1_001.1673523141941_Cut_0.fastq.gz"
    ,"/home/ubuntu/data/local/proj/bioinformatics_project/data/interim/exome/trimmed_umis/PCB75-PDX_S5_L001_R2_001.1673523141941_Cut_0.fastq.gz"
    ]
]

process BWA_align_human {
    publishDir "${params.outdir}/sam", mode: 'copy'
    maxForks 3
    input:
        path EXOME_trimmed_read_pair
    output:
        path "*.bam"
    """
    filename_human=\$(basename ${EXOME_trimmed_read_pair[0]} _Cut_0.fastq.gz)
    echo ${EXOME_trimmed_read_pair}
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_human} ${EXOME_trimmed_read_pair[0]} ${EXOME_trimmed_read_pair[1]} > ./\${filename_human::-21}_human.sam
    samtools view -bS \${filename_human::-21}_human.sam > \${filename_human::-21}.bam
    rm -f \${filename_human::-21}_human.sam 

    """
}


process BWA_align_mouse {
    publishDir "${params.outdir}", mode: 'copy'
    maxForks 4
    input:
        path EXOME_trimmed_read_pair
    output:
        file "*.bam"

    //https://biology.stackexchange.com/questions/59493/how-to-convert-bwa-mem-output-to-bam-format-without-saving-sam-file
    """
    filename_mouse=\$(basename ${EXOME_trimmed_read_pair[0]} _Cut_0.fastq.gz)
    bwa mem -t 7 -K 100000000 -Y ${params.fasta_mouse} ${EXOME_trimmed_read_pair[0]} ${EXOME_trimmed_read_pair[1]} | samtools sort -@7 -n -O BAM -o \${filename_mouse::-22}_mouse.sorted_by_name.bam -
    
    """
}


process SAM_sort_name {
    maxForks 9
    publishDir "${params.outdir}/sam_sorted", mode: 'copy'
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
    //trimmed_read_pairs_ch=Channel.fromFilePairs(params.EXOME_trimmed_reads_directory, flat:true)
    //trimmed_read_pairs_ch.view();
    //DISAMBIGUATE_human=BWA_align_human(trimmed_read_pairs_ch)
    //DISAMBIGUATE_mouse=BWA_align_mouse(trimmed_read_pairs_ch)
    //DISAMBIGUATE_ch=DISAMBIGUATE_human.merge(DISAMBIGUATE_mouse)
    //DISAMBIGUATE_ch.view()

    
   
   
    DISAMBIGUATE_all=Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/sam_sorted/*.bam")
    DISAMBIGUATE_ch=DISAMBIGUATE_all.map{it ->[it.name.split('_')[0],it] }.groupTuple()
   
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

