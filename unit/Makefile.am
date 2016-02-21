if BUNIT_TEST_LIB_PRESENT

check_PROGRAMS=mutateTest \
	crossoverTest \
	gameteTest \
	mutationTest \
	mlocusCrossoverTest \
	siteDepFitness \
	sugar_singlepop \
	sugar_metapop \
	sugar_multilocus \
	sugar_singlepop_custom_diploid \
	ms_sampling \
	sugar_sampling \
	callbacks \
	sugar_metapop_custom_diploid \
	test_generalmut \
	test_generalmut_vec \
	extensions \
	demography \
	gamete_cleaner \
	GSLrng_test \
	type_traits_test


TESTS=$(check_PROGRAMS)

mutateTest_SOURCES=mutateTest.cc
crossoverTest_SOURCES=crossoverTest.cc
gameteTest_SOURCES=gameteTest.cc
mutationTest_SOURCES=mutationTest.cc
mlocusCrossoverTest_SOURCES=mlocusCrossoverTest.cc
siteDepFitness_SOURCES=siteDepFitness.cc
sugar_singlepop_SOURCES=sugar_singlepop.cc
sugar_metapop_SOURCES=sugar_metapop.cc
sugar_multilocus_SOURCES=sugar_multilocus.cc
sugar_singlepop_custom_diploid_SOURCES=sugar_singlepop_custom_diploid.cc
sugar_metapop_custom_diploid_SOURCES=sugar_metapop_custom_diploid.cc
ms_sampling_SOURCES=ms_sampling.cc
sugar_sampling_SOURCES=sugar_sampling.cc
callbacks_SOURCES=callbacks.cc
test_generalmut_SOURCES=test_generalmut.cc
test_generalmut_vec_SOURCES=test_generalmut_vec.cc
extensions_SOURCES=extensions.cc
demography_SOURCES=demography.cc
gamete_cleaner_SOURCES=gamete_cleaner.cc
GSLrng_test_SOURCES=GSLrng_test.cc
type_traits_test_SOURCES=type_traits_test.cc

AM_CXXFLAGS=

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

endif
