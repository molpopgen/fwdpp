if BUNIT_TEST_PRESENT

noinst_PROGRAMS=unit/fwdpp_unit_tests unit/extensions_unit_tests unit/sugar_unit_tests \
				unit/genetic_map_tests \
				integration/sugar_integration_tests \
				tree_sequences/tree_sequence_tests

# integration/extensions_integration_tests 

TESTS=$(noinst_PROGRAMS)

#Unit test targets:
unit_fwdpp_unit_tests_SOURCES=unit/fwdpp_unit_tests.cc unit/mutateTest.cc unit/gameteTest.cc unit/utilTest.cc unit/type_traitsTest.cc unit/siteDepFitnessTest.cc unit/serializationTest.cc unit/gamete_cleanerTest.cc \
	unit/test_mutate_recombine.cc unit/test_sampling.cc \
	unit/test_simparams.cc \
	unit/test_wrapped_range.cc \
	unit/test_recombination.cc \
	fixtures/sugar_fixtures.cc

unit_genetic_map_tests_SOURCES=unit/genetic_map_tests.cc \
	unit/test_genetic_map.cc

unit_extensions_unit_tests_SOURCES=unit/extensions_unit_test.cc unit/extensions_regionsTest.cc unit/extensions_callbacksTest.cc
unit_sugar_unit_tests_SOURCES=unit/sugar_unit_tests.cc \
	unit/sugar_GSLrngTest.cc \
	unit/sugar_add_mutationTest.cc \
	unit/sugar_change_neutralTest.cc \
	unit/sugar_construct_pops.cc \
	unit/sugar_popgenmut.cc \
	fixtures/sugar_fixtures.cc

#Integration test targets:
integration_sugar_integration_tests_SOURCES=integration/sugar_integration_tests.cc integration/sugar_singlepop_custom_diploidTest.cc  integration/sugar_singlepopTest.cc \
	fixtures/sugar_fixtures.cc

# TODO: the next two test suites need to be brought back.

# integration_extensions_integration_tests_SOURCES=integration/extensions_integration_tests.cc integration/extensions_regionsIntegrationTest.cc \
# 	fixtures/sugar_fixtures.cc

tree_sequences_tree_sequence_tests_SOURCES=tree_sequences/tree_sequence_tests.cc \
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
										tree_sequences/test_ancestry_list.cc

AM_CXXFLAGS=-W -Wall

#AM_LIBS=-lboost_unit_test_framework

#LIBS+=$(AM_LIBS)

endif

