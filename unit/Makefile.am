if BUNIT_TEST_PRESENT

check_PROGRAMS=mutateTest \
	crossoverTest \
	gameteTest \
	mutationTest \
	policyTests \
	mlocusCrossoverTest

mutateTest_SOURCES=mutateTest.cc 
crossoverTest_SOURCES=crossoverTest.cc
gameteTest_SOURCES=gameteTest.cc
mutationTest_SOURCES=mutationTest.cc
policyTests_SOURCES=policyTests.cc
mlocusCrossoverTest_SOURCES=mlocusCrossoverTest.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif
