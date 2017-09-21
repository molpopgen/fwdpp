#!/usr/bin/env bash

if [ "$USECONDA" == "1" ];
then
    export LD_LIBRARY_PATH=$HOME/miniconda/lib
    export CPPFLAGS="-I$HOME/miniconda/include $CPPFLAGS"
    export LDFLAGS="-L$HOME/miniconda/lib $LDFLAGS"
    LDFLAGS="-L$HOME/miniconda/lib -Wl,-rpath,$HOME/miniconda/lib" ./configure --prefix=$HOME && make -j 3 &&  make install
    else
    CXXFLAGS="-std=c++11 -O2"
    if [ "$USECLANG" == "1" ];
    then
        CXXFLAGS+=" -stdlib=libc++"
    fi
    ./configure CXXFLAGS="$CXXFLAGS" --prefix=$HOME && make -j 3 && make install
fi

fwdppConfig --version
make check
