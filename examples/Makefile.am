noinst_PROGRAMS=	diploid_ind \
	diploid_fixed_sh_ind \
	K_linked_regions_extensions \
	juvenile_migration \
	custom_diploid_example \
	load_table_collection \
	edge_buffering \
	genetic_map_example \
	rng_example \
	diploid_example

diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
K_linked_regions_extensions_SOURCES=K_linked_regions_extensions.cc common_ind.hpp
juvenile_migration_SOURCES=juvenile_migration.cc common_ind.hpp

# These examples are used to test for API regressions and to include
# in the manual.
genetic_map_example_SOURCES=genetic_map_example.cc
rng_example_SOURCES=rng_example.cc
diploid_example_SOURCES=diploid_example.cc
custom_diploid_example_SOURCES=custom_diploid_example.cc

load_table_collection_SOURCES=load_table_collection.cc
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


edge_buffering_SOURCES=edge_buffering.cc \
					   $(tskit_SOURCES)

edge_buffering_CFLAGS=-std=c99
# examples based on tree sequences that need boost program options
if BPO_PRESENT
if BPO_PRESENT
noinst_PROGRAMS+=wfts_integration_test
wfts_integration_test_SOURCES=wfts_integration_test.cc tree_sequence_examples_common.cc
wfts_integration_test_LDADD=-lboost_program_options
endif
endif

edge_buffering_LDADD=-lboost_program_options

AM_CPPFLAGS=-Wall -W -I. -I../subprojects/nongpl/tskit/c -I../subprojects/nongpl/tskit/c/subprojects/kastore

AM_CXXFLAGS=
if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=
LIBS+=$(AM_LIBS)
