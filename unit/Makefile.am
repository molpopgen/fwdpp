if BUNIT_TEST_LIB_PRESENT

check_PROGRAMS=mutateTest \
	crossoverTest \
	gameteTest \
	mutationTest \
	mlocusCrossoverTest \
	siteDepFitness \
	sugar1 \
	sugar2 \
	sugar3 \
	sugar1_custom \
	ms_sampling \
	sugar_sampling \
	callbacks \
	sugar2_custom \
	test_generalmut \
	test_generalmut_vec


TESTS=$(check_PROGRAMS)

mutateTest_SOURCES=mutateTest.cc 
crossoverTest_SOURCES=crossoverTest.cc
gameteTest_SOURCES=gameteTest.cc
mutationTest_SOURCES=mutationTest.cc
policyTests_SOURCES=policyTests.cc
mlocusCrossoverTest_SOURCES=mlocusCrossoverTest.cc
siteDepFitness_SOURCES=siteDepFitness.cc
sugar1_SOURCES=sugar1.cc
sugar2_SOURCES=sugar2.cc
sugar3_SOURCES=sugar3.cc
sugar1_custom_SOURCES=sugar1_custom.cc
sugar2_custom_SOURCES=sugar2_custom.cc
ms_sampling_SOURCES=ms_sampling.cc
sugar_sampling_SOURCES=sugar_sampling.cc
callbacks_SOURCES=callbacks.cc
test_generalmut_SOURCES=test_generalmut.cc
test_generalmut_vec_SOURCES=test_generalmut_vec.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif
