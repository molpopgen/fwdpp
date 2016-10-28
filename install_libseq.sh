#!/usr/bin/env sh

##Install libseq 1.9.0
wget http://github.com/molpopgen/libsequence/archive/1.9.0.tar.gz
tar xzf 1.9.0.tar.gz
cd libsequence-1.9.0
./configure --prefix=$HOME/miniconda && make install
