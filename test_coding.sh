##run script

cd ~/nextflow_scripting/playground/test_pipeline

nextflow run  ~/00_scripts/nextflow/nucMACC --test -resume

nextflow run  ~/00_scripts/nextflow/nucMACC --blacklist ${baseDir}/toyData/dm3-blacklist_Chromosomes.bed -resume

nextflow run  ~/00_scripts/nextflow/nucMACC --test --analysis MNaseQC


