##run script

cd ~/nextflow_scripting/playground/test_pipeline

nextflow run  ~/00_scripts/nextflow/nucMACC --test -resume

nextflow run  ~/00_scripts/nextflow/nucMACC --blacklist ${baseDir}/toyData/dm3-blacklist_Chromosomes.bed -resume


### build docker

cd ~/00_scripts/nextflow/nucMACC/docker/common

docker build -t uschwartz/nucmacc:latest .


cd ~/00_scripts/nextflow/nucMACC/docker/qualimap

docker build -t uschwartz/qualimap:latest .



cd ~/00_scripts/nextflow/nucMACC/docker/danpos

docker build -t uschwartz/danpos:latest .



cd ~/00_scripts/nextflow/nucMACC/docker/R_nucMACC

docker build -t uschwartz/r_nucmacc:latest .


docker run --rm -it --volume /Users/admin/00_scripts/nextflow/nucMACC/docker/R_nucMACC:/Users/admin uschwartz/r_nucmacc:latest Rscript /Users/admin/hello.R

docker run --rm -it  uschwartz/r_nucmacc:latest
