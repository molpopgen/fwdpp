if HAVE_SIMDATA_HPP
if HAVE_LIBSEQ_RUNTIME

noinst_PROGRAMS=	diploid_ind \
	diploid_binaryIO_ind \
	diploid_gzbinaryIO_ind \
	diploid_fixed_sh_ind \
	migsel_ind \
	migsel_split_ind_list \
	bneck_selection_ind \
	diploid_ind_2locus \
	K_linked_regions_multilocus \
	K_linked_regions_extensions \
	HOC_ind 

#	migsel_split_ind 

diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_binaryIO_ind_SOURCES=diploid_binaryIO_ind.cc common_ind.hpp
diploid_gzbinaryIO_ind_SOURCES=diploid_gzbinaryIO_ind.cc common_ind.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
migsel_ind_SOURCES=migsel_ind.cc common_ind.hpp
migsel_split_ind_list_SOURCES=migsel_split_ind_list.cc common_ind.hpp
bneck_selection_ind_SOURCES=bneck_selection_ind.cc common_ind.hpp
diploid_ind_2locus_SOURCES=diploid_ind_2locus.cc common_ind.hpp
K_linked_regions_multilocus_SOURCES=K_linked_regions_multilocus.cc common_ind.hpp
K_linked_regions_extensions_SOURCES=K_linked_regions_extensions.cc common_ind.hpp
HOC_ind_SOURCES=HOC_ind.cc

AM_CPPFLAGS=-Wall -W -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=
AM_LIBS=-lsequence
LIBS+=$(AM_LIBS)
endif
endif
