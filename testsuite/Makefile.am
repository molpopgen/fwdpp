if BUNIT_TEST_PRESENT

check_PROGRAMS=unit/fwdpp_unit_tests unit/extensions_unit_tests unit/sugar_unit_tests \
			   integration/sugar_integration_tests \
			   integration/extensions_integration_tests

TESTS=$(check_PROGRAMS)

#Unit test targets:
unit_fwdpp_unit_tests_SOURCES=unit/fwdpp_unit_tests.cc unit/mutateTest.cc unit/gameteTest.cc unit/utilTest.cc unit/type_traitsTest.cc unit/siteDepFitnessTest.cc unit/serializationTest.cc unit/mlocusCrossoverTest.cc unit/gamete_cleanerTest.cc unit/test_general_rec_variation.cc \
	unit/test_mutate_recombine.cc unit/test_sampling.cc

unit_extensions_unit_tests_SOURCES=unit/extensions_unit_test.cc unit/extensions_regionsTest.cc unit/extensions_callbacksTest.cc
unit_sugar_unit_tests_SOURCES=unit/sugar_unit_tests.cc \
	unit/sugar_GSLrngTest.cc \
	unit/sugar_add_mutationTest.cc \
	unit/sugar_change_neutralTest.cc \
	unit/sugar_construct_pops.cc \
	unit/sugar_popgenmut.cc

#Integration test targets:
integration_sugar_integration_tests_SOURCES=integration/sugar_integration_tests.cc integration/sugar_multilocusTest.cc  integration/sugar_singlepop_custom_diploidTest.cc  integration/sugar_singlepopTest.cc

integration_extensions_integration_tests_SOURCES=integration/extensions_integration_tests.cc integration/extensions_regionsIntegrationTest.cc

AM_CXXFLAGS=-W -Wall

#AM_LIBS=-lboost_unit_test_framework

#LIBS+=$(AM_LIBS)

endif

