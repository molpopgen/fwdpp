if BUNIT_TEST_LIB_PRESENT

check_PROGRAMS=sugar_integration_tests


TESTS=$(check_PROGRAMS)

sugar_integration_tests_SOURCES=sugar_integration_tests.cc  sugar_metapop_custom_diploidTest.cc  sugar_metapopTest.cc  sugar_multilocusTest.cc  sugar_singlepop_custom_diploidTest.cc  sugar_singlepopTest.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif

