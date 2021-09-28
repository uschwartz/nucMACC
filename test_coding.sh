##run script

cd ~/nextflow_scripting/playground/test_pipeline

nextflow run  ~/00_scripts/nextflow/nucMACC -resume

nextflow run  ~/00_scripts/nextflow/nucMACC --blacklist ${baseDir}/toyData/dm3-blacklist_Chromosomes.bed -resume


### build docker

cd ~/00_scripts/nextflow/nucMACC/docker/common

docker build -t uschwartz/nucmacc:latest .


cd ~/00_scripts/nextflow/nucMACC/docker/qualimap

docker build -t uschwartz/qualimap:latest .

cd ~/00_scripts/nextflow/nucMACC/docker/danpos

docker build -t uschwartz/danpos:latest .
