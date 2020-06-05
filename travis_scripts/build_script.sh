#!/usr/bin/env bash

if [ "$USECONDA" == "1" ];
then
    echo "details..."
    echo `g++ -v`
    echo `gcc -v`
    echo "details done..."
    export PATH="$HOME/miniconda/bin:$PATH"
    if [ "$CXXSTANDARD" == "c++17" ];
    then
    CXX=g++ CC=gcc CPPFLAGS=-I$HOME/miniconda/include LDFLAGS="-L$HOME/miniconda/lib -Wl,-rpath,$HOME/miniconda/lib" ./configure --enable-cpp17 --prefix=$HOME && make -j 3 &&  make install
    else
    CXX=g++ CC=gcc CPPFLAGS=-I$HOME/miniconda/include LDFLAGS="-L$HOME/miniconda/lib -Wl,-rpath,$HOME/miniconda/lib" ./configure --prefix=$HOME && make -j 3 &&  make install
    fi
else
    if [ "$CXXSTANDARD" == "c++17" ];
    then
        ./configure --enable-cpp17 CXXFLAGS="$CXXFLAGS" --prefix=$HOME && make -j 3 && make install
    else
        ./configure CXXFLAGS="$CXXFLAGS" --prefix=$HOME && make -j 3 && make install
    fi
fi

fwdppConfig --version
make check
