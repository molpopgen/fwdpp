#!/usr/bin/env bash

echo "Standard = "

if [ "$USECONDA" == "1" ];
then
    export LD_LIBRARY_PATH=$HOME/miniconda/lib
    export CPPFLAGS="-I$HOME/miniconda/include $CPPFLAGS"
    export LDFLAGS="-L$HOME/miniconda/lib $LDFLAGS"
    LDFLAGS="-L$HOME/miniconda/lib -Wl,-rpath,$HOME/miniconda/lib" ./configure --prefix=$HOME && make -j 3 &&  make install
    else
    CXXFLAGS="-std=$CXXSTANDARD -O2"
    ./configure CXXFLAGS="$CXXFLAGS" --prefix=$HOME && make -j 3 && make install
fi

fwdppConfig --version
make check
./testsuite/unit/fwdpp_unit_tests --log_level=all
