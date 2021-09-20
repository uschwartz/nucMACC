##run script

cd ~/nextflow_scripting/playground/test_pipeline

nextflow run  ~/00_scripts/nextflow/nucMACC


### build docker

cd ~/00_scripts/nextflow/nucMACC/

docker build -t uschwartz/nucmacc:latest .
