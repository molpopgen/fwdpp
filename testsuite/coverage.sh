#!/bin/bash

gcov unit/*.cc integration/*.cc tree_sequences/*.cc
lcov -c -d . -o coverage_unfiltered.info
lcov -e coverage_unfiltered.info -o coverage_remove_boost.info '*fwdpp*'
lcov -r coverage_remove_boost.info -o coverage.info '*fwdpp/testsuite/*'
genhtml coverage.info -o ../fwdpp_coverage
