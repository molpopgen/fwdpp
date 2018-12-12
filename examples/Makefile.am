noinst_PROGRAMS=	diploid_ind \
	diploid_fixed_sh_ind \
	K_linked_regions_multilocus \
	K_linked_regions_extensions \
	K_linked_regions_generalized_rec \
	juvenile_migration \
	custom_mutation \
	custom_diploid


diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
K_linked_regions_multilocus_SOURCES=K_linked_regions_multilocus.cc common_ind.hpp
K_linked_regions_extensions_SOURCES=K_linked_regions_extensions.cc common_ind.hpp
K_linked_regions_generalized_rec_SOURCES=K_linked_regions_generalized_rec.cc common_ind.hpp
juvenile_migration_SOURCES=juvenile_migration.cc common_ind.hpp
custom_mutation_SOURCES=custom_mutation.cc common_ind.hpp
custom_diploid_SOURCES=custom_diploid.cc common_ind.hpp

# examples based on tree sequences
if BPO_PRESENT
if BPO_PRESENT
noinst_PROGRAMS+=wfts
wfts_SOURCES=wfts.cc tree_sequence_examples_common.cc
wfts_LDADD=-lboost_program_options

noinst_PROGRAMS+=wfmlts
wfmlts_SOURCES=wfmlts.cc tree_sequence_examples_common.cc
wfmlts_LDADD=-lboost_program_options

noinst_PROGRAMS+=spatialts
spatialts_SOURCES=spatialts.cc
spatialts_LDADD=-lboost_program_options
endif
endif

AM_CPPFLAGS=-Wall -W -I.

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
