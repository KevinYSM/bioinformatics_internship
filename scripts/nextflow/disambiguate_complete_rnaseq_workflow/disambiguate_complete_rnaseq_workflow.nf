#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow"


params.RNA_aligned_directory="/data/local/proj/bioinformatics_project/data/interim/rna/aligned"
params.RNA_raw_reads_directory="/data/local/proj/bioinformatics_project/data/raw/rna/PCB-28-PDX_S4_L004_R{1,2}_*.fastq.gz"
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
        memory '40G'
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
        memory '40G'
        input:
                file RNA_read_pair


        output:
        path '*.bam'

        """
        STAR_align.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]} ${params.fasta_mouse} ${params.gtf_mouse} ${params.mouseGenomeDir} 'mouse'
        """
}

process SAM_sort_name{
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

process SAM_sort{
     
        input:
        file DISAMBIGUATE_bam_file

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
        publishDir params.outdir, mode: 'copy'
        input:
        file SORTED_bam_file

        output:
        path '*.gene_counts'

        """
        HTSEQ_count.sh $SORTED_bam_file
        """
}
process GATK_mark_duplicates{
        input:
        file SORTED_bam_file

        output:
        path '*.bam'

        """
        GATK_mark_duplicates.sh $SORTED_bam_file 
        """
}

process GATK_split{
        input:
        file DUPLICATES_bam_file

        output:
        path '*.bam'
        """
        GATK_split_N_cigar_reads.sh $DUPLICATES_bam_file ${params.fasta_human}
        """
}

process GATK_base_recalibration{
        publishDir params.outdir, mode:'copy'
        input:
        file SPLIT_bam_file

        output:
        path '*'

        """
        GATK_base_recalibrator.sh $SPLIT_bam_file ${params.fasta_human} ${params.k1} ${params.k2}
        """
}


workflow{
        read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        STAR_ch_human=STAR_human(read_pairs_ch)
        STAR_ch_mouse=STAR_mouse(read_pairs_ch)

        DISAMBIGUATE_human=SAM_sort_name(STAR_ch_human)
        DISAMBIGUATE_mouse=SAM_sort_name_2(STAR_ch_mouse)

        DISAMBIGUATE_ch=DISAMBIGUATE_human.merge(DISAMBIGUATE_mouse)
        SAM_sort_ch=NGS_disambiguate(DISAMBIGUATE_ch)
        GATK_duplicates_ch=SAM_sort(SAM_sort_ch)

        HTSEQ_count(SAM_sort_ch)

        GATK_split_ch=GATK_mark_duplicates(GATK_duplicates_ch)
        GATK_split_ch.view()
        GATK_base_recalibration_ch=GATK_split(GATK_split_ch)
        GATK_base_recalibration_ch.view()
        GATK_base_recalibration(GATK_base_recalibration_ch)
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