#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/"


//params.RNA_aligned_directory="/data/local/proj/bioinformatics_project/data/interim/rna/aligned"
params.RNA_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/rna/*_L004_R{1,2}_*.fastq.gz"
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna"
params.gtf_human="/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
//params.gtf_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.gtf.gz"
params.gtf_mouse="/data/local/reference/mouse/test/GCF_000001635.27_GRCm39_genomic.gtf"
params.humanGenomeDir="/data/local/proj/bioinformatics_project/data/external/human_genome"
params.mouseGenomeDir="/data/local/proj/bioinformatics_project/data/external/mouse_genome"

params.k1="/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz"
params.k2="/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

process STAR_human{
        maxForks 1

        input:
                file RNA_read_pair
        output:
        path '*.bam'

        """
        STAR_align.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]} ${params.fasta_human} ${params.gtf_human} ${params.humanGenomeDir} 'human'
        """
}

process STAR_mouse{
        maxForks 1
        input:
                file RNA_read_pair


        output:
        path '*.bam'

        """
        STAR_align.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]} ${params.fasta_mouse} ${params.gtf_mouse} ${params.mouseGenomeDir} 'mouse'
        """
}

process SAM_sort_name{
        publishDir "${params.outdir}/SAM_sort_name" , mode: 'symlink'
        input:
        file ALIGNED_bam_file

        output:
        path '*.bam'

        """
        SAM_sort_name.sh $ALIGNED_bam_file
        """
}

process SAM_sort_name_2{
        input:
        file ALIGNED_bam_file

        output:
        path '*.bam'

        """
        SAM_sort_name.sh $ALIGNED_bam_file
        """
}



process NGS_disambiguate{
        publishDir "${params.outdir}/disambiguate/" , mode: 'copy'

        input:
                tuple val(key), path(BAM_files)

        output:
                file '*.disambiguatedSpeciesA.bam'
                file '*.txt'

        """
       if [ ${BAM_files[0]} == *"human"* ]
        then
                echo "human first"
                echo ${BAM_files}
                NGS_disambiguate.sh ${BAM_files[0]} ${BAM_files[1]}
               
        elif [ ${BAM_files[0]} == *"mouse"* ]      
        then
                echo "mouse first"
                echo ${BAM_files}
                NGS_disambiguate.sh ${BAM_files[1]} ${BAM_files[0]}

        else
                echo "Error"
        fi

        """
}

process SAM_sort{
        publishDir "${params.outdir}/SAM_sort" , mode: 'symlink'
        input:
                file DISAMBIGUATE_bam_file
                file DISAMBIGUATED_txt

        output:
                path '*.bam'


        """
        SAM_sort.sh $DISAMBIGUATE_bam_file

        
        """
}


process SAM_sort_name_3{
        input:
        file ALIGNED_bam_file

        output:
                path '*.bam'

        """
        SAM_sort_name.sh $ALIGNED_bam_file
        """

}




process HTSEQ_count{
        maxForks 3
                publishDir "${params.outdir}/gene_counts", mode: 'copy'
        input:
                file SORTED_bam_file

        output:
                file '*.gene_counts'

        """
        HTSEQ_count.sh $SORTED_bam_file
        """
}
process GATK_mark_duplicates{
         maxForks 3
        publishDir "${params.outdir}/mark_duplicates", mode: 'symlink'
        input:
                file SORTED_bam_file

        output:
                file '*.bam'

        """
        GATK_mark_duplicates.sh $SORTED_bam_file 
        """
}

process GATK_split{
      maxForks 3
        publishDir "${params.outdir}/split", mode: 'symlink'
        
        input:
                file DUPLICATES_bam_file

        output:
        path '*.bam'
        """
        GATK_split_N_cigar_reads.sh $DUPLICATES_bam_file ${params.fasta_human}
        """
}

process GATK_base_recal_all{
        maxForks 3
        publishDir "${params.outdir}", mode: 'copy'
     
        input:
                file SPLIT_bam_file
        output:
                file "*.recal.pass2.bam"
                file "*.recal.pass2.table"
                file "*.pass2.bai"
        """
        GATK_base_recal_all.sh $SPLIT_bam_file ${params.fasta_human} ${params.k1} ${params.k2}
        """
}

process GATK_base_recalibrator_1{
        maxForks 2
        input:
                file SPLIT_bam_file

        output:
                path '*.recal.pass1.table'
                

        """
        GATK_base_recalibrator_1.sh $SPLIT_bam_file ${params.fasta_human} ${params.k1} ${params.k2}
        """
}
process GATK_apply_bqsr_1 {
        maxForks 2
        publishDir "${params.outdir}" , mode: 'symlink'
        input:
        file BASE_recalibrated_table

        output:
        path '*.recal.pass1.bam'

        """
        GATK_apply_bqsr_1.sh $BASE_recalibrated_table ${params.fasta_human} ${params.k1} ${params.k2}
        """
}

process GATK_base_recal_and_bqsr_1{
        input:
                file SPLIT_bam_file

        output:
                path '*.recal.pass1.bam'

        """
        GATK_base_recalibrator_1.sh $SPLIT_bam_file ${params.fasta_human} ${params.k1} ${params.k2} | GATK_apply_bqsr_1.sh 
        """
}

process GATK_base_recalibrator_2 {
        maxForks 2
        input:
        file BAM_bqsr_file

        output:
        path '*.recal.pass2.table'

        """
        GATK_base_recalibrator_2.sh $BAM_bqsr_file ${params.fasta_human} ${params.k1} ${params.k2}
        """
}

process GATK_apply_bqsr_2 {
        maxForks 2
        publishDir "${params.outdir}", mode: 'copy'
        input:
        file BASE_recalibrated_file

        output:
        path '*'

        """
        GATK_apply_bqsr_2.sh $BASE_recalibrated_file ${params.fasta_human} ${params.k1} ${params.k2}
        """

        
}
workflow{
        read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)

        //Perform Alignment
        ///STAR_ch_human=STAR_human(read_pairs_ch)
        ///STAR_ch_mouse=STAR_mouse(read_pairs_ch)

        //Sort all bam files
        ///SAM_sort_name_ch=STAR_ch_human.mix(STAR_ch_mouse)

        //Sort bam files by name
        ///SORTED_bam_files_ch=SAM_sort_name(SAM_sort_name_ch)

        //Create .bam human and mouse tuples
        ///SORTED_bam_files_ch=Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/SAM_sort_name/*.bam" )
        ///DISAMBIGUATE_ch=SORTED_bam_files_ch.map{it ->[it.name.split('_')[0],it] }.groupTuple()

        // have to add .collect before .map
        ///DISAMBIGUATE_ch.view()

        //Perform Disambiguate on channel
        ///SAM_sort_ch=NGS_disambiguate(DISAMBIGUATE_ch)

        //Sort disambiguated bam file by index
        ///GATK_duplicates_ch=SAM_sort(SAM_sort_ch)
        ///GATK_duplicates_ch=Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/SAM_sort/PCB-39-PDX_S8_L004.disambiguatedSpeciesA.sorted.bam")
        //a) Perform HTSeq using Sorted Channel
        ///HTSEQ_count(GATK_duplicates_ch)

        //Mark duplicate reads
        ///GATK_split_ch=GATK_mark_duplicates()
 
        //Split N Cigar Reads
        ///GATK_base_recalibration_ch=GATK_split(Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/mark_duplicates/*.bam"))

        //Perform Two Rounds of Base Recalibration
        GATK_base_recal_all(Channel.fromPath("/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow/split/PCB-09*.split.bam"))


        //HTSEQ_count_ch=DISAMBIGUATE_human.concat(DISAMBIGUATE_mouse)

        //GATK_duplicates_ch=SAM_index(SAM_index_ch)
        
        //read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        //STAR_ch=STAR_align_human(read_pairs_ch)
        //SAM_ch=SAM_sort_name(STAR_ch)
        //HTSEQ_count(SAM_ch)
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