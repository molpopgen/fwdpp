if BUNIT_TEST_PRESENT 
if BUNIT_TEST_RUNTIME_PRESENT

noinst_PROGRAMS=unit/fwdpp_unit_tests unit/extensions_unit_tests unit/sugar_unit_tests \
				unit/genetic_map_tests \
				integration/sugar_integration_tests \
				tree_sequences/tree_sequence_tests

# integration/extensions_integration_tests 

TESTS=$(noinst_PROGRAMS)

#Unit test targets:
unit_fwdpp_unit_tests_SOURCES=unit/fwdpp_unit_tests.cc unit/mutateTest.cc unit/gameteTest.cc unit/utilTest.cc unit/test_type_traits.cc unit/siteDepFitnessTest.cc unit/serializationTest.cc unit/gamete_cleanerTest.cc \
	unit/test_mutate_recombine.cc \
	unit/test_simparams.cc \
	unit/test_wrapped_range.cc \
	unit/test_enum_bitflags.cc \
	unit/test_nested_forward_lists.cc \
	unit/test_validators.cc \
	fixtures/sugar_fixtures.cc \
	fixtures/fwdpp_fixtures.hpp \
	fixtures/sugar_fixtures.hpp \
	util/custom_dip.hpp

unit_genetic_map_tests_SOURCES=unit/genetic_map_tests.cc \
	unit/test_genetic_map.cc \
	fixtures/rng_fixture.hpp

unit_extensions_unit_tests_SOURCES=unit/extensions_unit_test.cc unit/extensions_regionsTest.cc unit/extensions_callbacksTest.cc
unit_sugar_unit_tests_SOURCES=unit/sugar_unit_tests.cc \
	unit/sugar_GSLrngTest.cc \
	unit/sugar_add_mutationTest.cc \
	unit/sugar_change_neutralTest.cc \
	unit/sugar_popgenmut.cc \
	fixtures/sugar_fixtures.cc \
	fixtures/sugar_fixtures.hpp \
	util/custom_dip.hpp

#Integration test targets:
integration_sugar_integration_tests_SOURCES=integration/sugar_integration_tests.cc integration/sugar_singlepop_custom_diploidTest.cc  integration/sugar_singlepopTest.cc \
	fixtures/sugar_fixtures.cc \
	util/quick_evolve_sugar.hpp

# TODO: the next two test suites need to be brought back.

# integration_extensions_integration_tests_SOURCES=integration/extensions_integration_tests.cc integration/extensions_regionsIntegrationTest.cc \
# 	fixtures/sugar_fixtures.cc
tskit_SOURCES=../subprojects/nongpl/tskit/c/tskit/convert.c \
			  ../subprojects/nongpl/tskit/c/tskit/core.c \
			  ../subprojects/nongpl/tskit/c/tskit/genotypes.c \
			  ../subprojects/nongpl/tskit/c/tskit/haplotype_matching.c \
			  ../subprojects/nongpl/tskit/c/tskit/stats.c \
			  ../subprojects/nongpl/tskit/c/tskit/tables.c \
			  ../subprojects/nongpl/tskit/c/tskit/trees.c \
			  ../subprojects/nongpl/tskit/c/subprojects/kastore/kastore.c \
			  ../subprojects/nongpl/tskit/c/tskit.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/convert.h \
		  	  ../subprojects/nongpl/tskit/c/tskit/core.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/genotypes.h \
		  	  ../subprojects/nongpl/tskit/c/tskit/haplotype_matching.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/stats.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/tables.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/trees.h \
	  		  ../subprojects/nongpl/tskit/c/tskit/trees.h \
	  		  ../subprojects/nongpl/tskit/c/subprojects/kastore/kastore.h


tree_sequences_tree_sequence_tests_SOURCES=tree_sequences/tree_sequence_tests.cc \
										${tskit_SOURCES} \
										tree_sequences/empty_table_collection.hpp \
										tree_sequences/independent_implementations.hpp \
										tree_sequences/preorder_adl.hpp \
										tree_sequences/simple_table_collection.hpp \
										tree_sequences/simple_table_collection_infinite_sites.hpp \
										tree_sequences/simple_table_collection_polytomy.hpp \
										tree_sequences/tskit_utils.hpp \
										tree_sequences/wfevolve_table_collection.hpp \
										tree_sequences/test_generate_offspring.cc \
									    tree_sequences/test_preorder_node_traversal.cc \
										tree_sequences/test_tree_visitor.cc \
										tree_sequences/test_decapitate.cc \
										tree_sequences/test_marginal_tree_statistics.cc \
										tree_sequences/test_node_traversal_order_adl.cc \
										tree_sequences/test_sample_traversal.cc \
										tree_sequences/test_root_traversal.cc \
										tree_sequences/test_node_children_traversal.cc \
										tree_sequences/independent_implementations.cc \
										tree_sequences/test_generate_data_matrix.cc \
										tree_sequences/test_table_collection.cc \
										tree_sequences/test_visit_sites.cc \
										tree_sequences/test_site_visitor.cc \
										tree_sequences/test_marginal_tree.cc \
										tree_sequences/test_ancestry_list.cc \
										tree_sequences/test_diploid_recording.cc \
										tree_sequences/wfevolve_table_collection_fxns.cc \
										tree_sequences/test_edge_buffering_std_table_collection.cc \
										tree_sequences/tskit_utils.cc

tree_sequences_tree_sequence_tests_CFLAGS=-std=c99

AM_CPPFLAGS=-I../subprojects/nongpl/tskit/c -I../subprojects/nongpl/tskit/c/subprojects/kastore
AM_CXXFLAGS=-W -Wall --coverage -DBOOST_TEST_DYN_LINK

AM_LIBS=-lboost_unit_test_framework

LIBS+=$(AM_LIBS)

coverage-local: check coverage.sh
	sh coverage.sh

endif
endif
