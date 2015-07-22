#!/usr/bin/env sh

##Install libseq 1.8.5
wget http://github.com/molpopgen/libsequence/archive/1.8.5.tar.gz
tar xzf 1.8.5.tar.gz
cd libsequence-1.8.5
./configure --prefix=$HOME && make install
