#!/usr/bin/env bash

if [ "$USECONDA" == "1" ];
then
    export LD_LIBRARY_PATH=$HOME/miniconda/lib
    export CPPFLAGS="-I$HOME/miniconda/include $CPPFLAGS"
    export LDFLAGS="-L$HOME/miniconda/lib $LDFLAGS"
    LDFLAGS="-L$HOME/miniconda/lib -Wl,-rpath,$HOME/miniconda/lib" ./configure --prefix=$HOME && make -j 3 &&  make install
    else
    ./configure --prefix=$HOME && make && make install
fi

fwdppConfig --version
make check
