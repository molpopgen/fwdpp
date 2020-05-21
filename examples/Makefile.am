noinst_PROGRAMS=	diploid_ind \
	diploid_fixed_sh_ind \
	K_linked_regions_extensions \
	juvenile_migration \
	custom_mutation \
	custom_diploid \
	wfts_overlapping_generations \
	wfts_overlapping_generations_dynamic_indexing \
	wfts_overlapping_generations_dynamic_indexing_with_mutation \
	load_table_collection \
	edge_buffering


diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
K_linked_regions_extensions_SOURCES=K_linked_regions_extensions.cc common_ind.hpp
juvenile_migration_SOURCES=juvenile_migration.cc common_ind.hpp
custom_mutation_SOURCES=custom_mutation.cc common_ind.hpp
custom_diploid_SOURCES=custom_diploid.cc common_ind.hpp
wfts_overlapping_generations_SOURCES=wfts_overlapping_generations.cc wfevolvets.cc wfevolvets.hpp
wfts_overlapping_generations_dynamic_indexing_SOURCES=wfts_overlapping_generations_dynamic_indexing.cc wfevolvets.cc wfevolvets.hpp
wfts_overlapping_generations_dynamic_indexing_with_mutation_SOURCES=wfts_overlapping_generations_dynamic_indexing_with_mutation.cc wfevolvets.cc wfevolvets.hpp
load_table_collection_SOURCES=load_table_collection.cc
tskit_SOURCES=../subprojects/nongpl/tskit/c/tskit/convert.c \
			  ../subprojects/nongpl/tskit/c/tskit/core.c \
			  ../subprojects/nongpl/tskit/c/tskit/genotypes.c \
			  ../subprojects/nongpl/tskit/c/tskit/haplotype_matching.c \
			  ../subprojects/nongpl/tskit/c/tskit/stats.c \
			  ../subprojects/nongpl/tskit/c/tskit/tables.c \
			  ../subprojects/nongpl/tskit/c/tskit/trees.c \
			  ../subprojects/nongpl/tskit/c/subprojects/kastore/kastore.c
edge_buffering_SOURCES=edge_buffering.cc \
					   $(tskit_SOURCES)

edge_buffering_CFLAGS=-std=c99
# examples based on tree sequences that need boost program options
if BPO_PRESENT
if BPO_PRESENT
noinst_PROGRAMS+=wfts_integration_test
wfts_integration_test_SOURCES=wfts_integration_test.cc tree_sequence_examples_common.cc
wfts_integration_test_LDADD=-lboost_program_options

noinst_PROGRAMS+=spatialts
spatialts_SOURCES=spatialts.cc
spatialts_LDADD=-lboost_program_options
endif
endif

edge_buffering_LDADD=-lboost_program_options

AM_CPPFLAGS=-Wall -W -I. -I../subprojects/nongpl/tskit/c -I../subprojects/nongpl/tskit/c/subprojects/kastore

AM_CXXFLAGS=
if HAVE_SIMDATA_HPP
if HAVE_LIBSEQ_RUNTIME
AM_CXXFLAGS+=-DHAVE_LIBSEQUENCE
AM_LIBS=-lsequence
endif
endif
if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=
LIBS+=$(AM_LIBS)
