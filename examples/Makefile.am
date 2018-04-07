noinst_PROGRAMS=	diploid_ind \
	diploid_fixed_sh_ind \
	K_linked_regions_multilocus \
	K_linked_regions_extensions \
	K_linked_regions_generalized_rec 


diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
K_linked_regions_multilocus_SOURCES=K_linked_regions_multilocus.cc common_ind.hpp
K_linked_regions_extensions_SOURCES=K_linked_regions_extensions.cc common_ind.hpp
K_linked_regions_generalized_rec_SOURCES=K_linked_regions_generalized_rec.cc common_ind.hpp

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
