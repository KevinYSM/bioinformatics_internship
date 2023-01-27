params.outdir="/data/local/proj/bioinformatics_project/data/processed/sarek_workflow/final/"
params.fastqdir="data/local/data/interim/exome/trimmed_umis/"




process{
        nextflow run nf-core/sarek --input samplesheet.csv --outdir "data/local/proj/bioinformatics_project/data/processed/sarek_workflow/final/" -profile singularity --genome GTK.GRCh38 --cpus 31 --tools Mutect2,HaploTypeCaller,ASCAT,ControlFREEC,MSIsensor --wes
}
