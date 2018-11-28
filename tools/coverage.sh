#!/bin/bash

if [ -d coverage ]
then 
    rm -rf coverage
fi

if [ -e lcov.info ]
then
    rm -f lcov.info
fi
./configure CXXFLAGS="-O2 -Wall --coverage"
make clean
make -j 6
lcov -c -i -d testsuite --base-directory ./fwdpp --no-external -o base_lcov.info
make check -j 6
lcov -c -d testsuite --base-directory ./fwdpp --no-external -o test_lcov.info
lcov -a base_lcov.info -a test_lcov.info -o lcov.info
lcov --remove lcov.info '*/testsuite/*' '*/fwdpp/unit/*' '*/fwdpp/util/*' '*/fwdpp/fixtures/*' -o lcov.info
genhtml -s --ignore-errors source -o coverage lcov.info
