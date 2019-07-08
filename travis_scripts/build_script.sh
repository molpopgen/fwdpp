#!/usr/bin/env bash

if [ "$USECONDA" == "1" ];
then
    echo "details..."
    echo `g++ -v`
    echo $CXX
    echo "details done..."
    ./configure --prefix=$HOME && make -j 3 &&  make install
else
    CXXFLAGS="-std=$CXXSTANDARD -O2"
    ./configure CXXFLAGS="$CXXFLAGS" --prefix=$HOME && make -j 3 && make install
fi

fwdppConfig --version
make check
