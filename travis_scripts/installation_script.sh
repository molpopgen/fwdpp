#!/usr/bin/env bash

## This is a complex script...

echo "Variables = "$USECONDA $TRAVIS_OS_NAME

if [ "$USECONDA" == "1" ];
then
    if [ "$TRAVIS_OS_NAME" == "linux" ];
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [ "$TRAVIS_OS_NAME" == "osx" ];
        then
            wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    # Useful for debugging any issues with conda
    conda info -a
    conda install gxx_linux-64 zlib gsl
    conda install -c conda-forge boost-cpp==1.70.0
    ls -lhrt $HOME/miniconda/lib | grep boost_program
    # conda install -c bioconda libsequence
else
    sudo apt-get update -qq
    sudo apt-get -f install libboost-test-dev libboost-program-options-dev libgsl*-dev
fi
