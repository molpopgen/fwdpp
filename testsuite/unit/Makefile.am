if BUNIT_TEST_LIB_PRESENT

check_PROGRAMS=fwdpp_unit_tests


TESTS=$(check_PROGRAMS)

fwdpp_unit_tests_SOURCES=fwdpp_unit_tests.cc mutateTest.cc gameteTest.cc utilTest.cc type_traitsTest.cc recombinationTest.cc
#fwdpp_internal_unit_tests_SOURCES=fwdpp_internal_unit_tests.cc crossoverTest.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif

