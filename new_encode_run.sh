#!/bin/bash

codedir=/home/timp/Code

##Using the pipeline at https://github.com/ENCODE-DCC/atac-seq-pipeline

##To install genome
if true; then

    cd ${codedir}/atac-seq-pipeline/installers/
    mkdir -p human
    cp -f install_genome_data.sh human
    docker run -v $(cd $(dirname human) && pwd -P)/$(basename human):/genome_data_tmp quay.io/encode-dcc/atac-seq-pipeline:v1 "cd /genome_data_tmp && bash install_genome_data.sh hg38 ."

fi


if false; then
    cd ${codedir}/atac-seq-pipeline
    ##downloaded cromwell34 from https://github.com/broadinstitute/cromwell/releases/tag/34, using that instead of 30.2
    java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run atac.wdl -i /home/timp/Code/fragillis/fragillis_atac_24hrbft2.json -o workflow_opts/docker.json
fi
