params.SAM_directory="/data/local/proj/bioinformatics_project/data/interim/exome/sam/PCB{64,75}*"
params.GENDER_info_txt="/data/local/proj/bioinformatics_project/scripts/nextflow/gender_determination_workflow/gender_information.txt"
process SEX_determine{
    input:
        path SAM_file
    output:

    """
    ./sex_samples.sh ${SAM_file} ${params.GENDER_info_txt}
    """
    }

workflow{
    SEX_determine(Channel.fromPath(params.SAM_directory, flat:true))
}