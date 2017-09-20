#!/usr/bin/env sh

## This is a complex script...

echo $USECONDA $TRAVIS_OS_NAME

if [ "$USECONDA" == "1" ]
then
    if [ "$TRAVIS_OS_NAME" == "linux" ]
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [ "$TRAVIS_OS_NAME" == "osx" ]
        then
            wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    Useful for debugging any issues with conda
    conda info -a
    conda install gcc zlib boost
    conda install -c asmeurer gsl
    conda install -c bioconda libsequence
fi
