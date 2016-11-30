if BUNIT_TEST_PRESENT

noinst_PROGRAMS=unit/fwdpp_unit_tests unit/extensions_unit_tests unit/sugar_unit_tests \
			   integration/sugar_integration_tests

TESTS=$(noinst_PROGRAMS)

#Unit test targets:
unit_fwdpp_unit_tests_SOURCES=unit/fwdpp_unit_tests.cc unit/mutateTest.cc unit/gameteTest.cc unit/utilTest.cc unit/type_traitsTest.cc unit/recombinationTest.cc unit/demographyTest.cc unit/siteDepFitnessTest.cc unit/serializationTest.cc unit/ms_samplingTest.cc unit/mlocusCrossoverTest.cc unit/gamete_cleanerTest.cc
unit_extensions_unit_tests_SOURCES=unit/extensions_unit_test.cc unit/extensions_regionsTest.cc unit/extensions_callbacksTest.cc
unit_sugar_unit_tests_SOURCES=unit/sugar_unit_tests.cc \
	unit/sugar_GSLrngTest.cc \
	unit/sugar_samplingTest.cc \
	unit/sugar_add_mutationTest.cc \
	unit/sugar_change_neutralTest.cc \
	unit/sugar_generalmutTest.cc \
	unit/sugar_generalmut_vecTest.cc

#Integration test targets:
integration_sugar_integration_tests_SOURCES=integration/sugar_integration_tests.cc  integration/sugar_metapop_custom_diploidTest.cc  integration/sugar_metapopTest.cc  integration/sugar_multilocusTest.cc  integration/sugar_singlepop_custom_diploidTest.cc  integration/sugar_singlepopTest.cc integration/sugar_matrixTest.cc

AM_CXXFLAGS=

#AM_LIBS=-lboost_unit_test_framework

#LIBS+=$(AM_LIBS)

endif

