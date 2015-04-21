if BUNIT_TEST_PRESENT

check_PROGRAMS=mutateTest \
	crossoverTest \
	gameteTest \
	mutationTest \
	policyTests \
	mlocusCrossoverTest \
	siteDepFitness \
	sugar1 \
	sugar2

mutateTest_SOURCES=mutateTest.cc 
crossoverTest_SOURCES=crossoverTest.cc
gameteTest_SOURCES=gameteTest.cc
mutationTest_SOURCES=mutationTest.cc
policyTests_SOURCES=policyTests.cc
mlocusCrossoverTest_SOURCES=mlocusCrossoverTest.cc
siteDepFitness_SOURCES=siteDepFitness.cc
sugar1_SOURCES=sugar1.cc
sugar2_SOURCES=sugar2.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif
