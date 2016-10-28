if HAVE_SIMDATA_HPP
if HAVE_LIBSEQ_RUNTIME

noinst_PROGRAMS=sex_limited_ind diploid_fixed_sh_ind

sex_limited_ind_SOURCES=sex_limited_ind.cc 
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc

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
