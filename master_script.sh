#!/bin/bash

if [ "$1" == "conda.py27" ]; then

    conda create --name py27.fragillis
    source activate py27.fragillis
    conda update --all --yes
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda install macs2 haystack_bio numpy
    source deactivate
    

fi

if [ "$1" == "conda" ]; then

    conda create --name fragillis
    source activate fragillis
    conda install -c bioconda deeptools
    source deactivate

fi
    
