#!/usr/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// params.PROJECT_DIRECTORY="/data/local/proj/bioinformatics_project/data/processed/complete_rnaseq_workflow_nextflow"
params.outdir="/data/local/proj/bioinformatics_project/data/processed/disambiguate_complete_rnaseq_workflow_nextflow"


params.RNA_aligned_directory='/data/local/proj/bioinformatics_project/data/interim/rna/aligned'
params.RNA_raw_reads_directory='/data/local/proj/bioinformatics_project/data/raw/rna/*_R{1,2}_*.fastq.gz'
params.fasta_human="/data/local/reference/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
params.fasta_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.fna.gz"
params.gtf_human="/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"
params.gtf_mouse="/data/local/reference/mouse/GCF_000001635.27_GRCm39_genomic.gtf.gz"
params.humanGenomeDir="/data/local/proj/bioinformatics_project/data/external/human_genome"
params.mouseGenomeDir="/data/local/proj/bioinformatics_project/data/external/mouse_genome"

params.k1="/data/local/reference/GATK_resource_bundle/hg38/hg38/dbsnp_146.hg38.vcf.gz"
params.k2="/data/local/reference/GATK_resource_bundle/hg38/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

process STAR_human{
        input:
        file RNA_read_pair
        output:
        path '*.bam'

        """
        STAR_align.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]} ${params.fasta_human} ${params.gtf_human} ${params.humanGenomeDir}
        """


}

process STAR_mouse{
        input:
        file RNA_read_pair
        output:
        path '*.bam'

        """
        STAR_align.sh ${RNA_read_pair[1]} ${RNA_read_pair[2]} ${params.fasta_mouse} ${params.gtf_mouse} ${params.mouseGenomeDir}
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
        file BAM_files

        output:
        path '*.bam'

        """
      
        NGS_disambiguate.sh ${BAM_files[1]} ${BAM_files[2]}
        """

}

process SAM_sort{
     
        input:
        file DISAMBIGUATE_bam_file

        output:
        path '*.host.bam'

        """
        SAM_sort.sh $DISAMBIGUATE_bam_file
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

process GATK_split_ch{
        input:
        file DUPLICATES_bam_file

        output:
        path '*.bam'
        """
        GATK_split_N_cigar_reads.sh $DUPLICATES_bam_file ${params.fasta}
        """
}

process GATK_base_recalibration{
        publishDir params.outdir, mode:'copy'
        input:
        file SPLIT_bam_file

        output:
        path '*'

        """
        GATK_base_recalibrator.sh $SPLIT_bam_file ${params.fasta} ${k1} ${k2}
        """
}



workflow{
        read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        STAR_ch_human=STAR_human(read_pairs_ch)
        STAR_ch_mouse=STAR_mouse(read_pairs_ch)

        DISAMBIGUATE_human=SAM_sort_name(STAR_ch_human)
        DISAMBIGUATE_mouse=SAM_sort_name_2(STAR_ch_mouse)

        DISAMBIGUATE_ch=DISAMBIGUATE_human.join(DISAMBIGUATE_mouse)

        
        SAM_sort_ch=NGS_disambiguate(DISAMBIGUATE_ch)
        GATK_duplicates_ch=SAM_sort(SAM_sort_ch)


        //GATK_duplicates_ch=SAM_index(SAM_index_ch)
        GATK_split_ch=GATK_mark_duplicates(GATK_duplicates_ch)
        GATK_base_recalibration(GATK_split_ch)
        //read_pairs_ch = Channel.fromFilePairs(params.RNA_raw_reads_directory, flat: true)
        //STAR_ch=STAR_align_human(read_pairs_ch)
        //SAM_ch=SAM_sort_name(STAR_ch)
        //HTSEQ_count(SAM_ch)
        
}