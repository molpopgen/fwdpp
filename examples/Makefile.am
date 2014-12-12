check_PROGRAMS=diploid \
	diploid_ind \
	diploid_binaryIO \
	diploid_binaryIO_ind \
	diploid_gzbinaryIO \
	diploid_gzbinaryIO_ind \
	diploid_fixed_sh \
	diploid_fixed_sh_ind \
	diploid_fixed_sh_ind_lambda \
	migsel_ind \
	bneck_selection_ind \
	diploid_ind_2locus \
	diploid_twopop_mig \
	pfix \
	RHH \
	TFL2013 

diploid_SOURCES=diploid.cc common_gamete.hpp
diploid_ind_SOURCES=diploid_ind.cc common_ind.hpp
diploid_binaryIO_SOURCES=diploid_binaryIO.cc common_gamete.hpp
diploid_binaryIO_ind_SOURCES=diploid_binaryIO_ind.cc common_ind.hpp
diploid_gzbinaryIO_SOURCES=diploid_gzbinaryIO.cc common_gamete.hpp
diploid_gzbinaryIO_ind_SOURCES=diploid_binaryIO_ind.cc common_ind.hpp
diploid_fixed_sh_SOURCES=diploid_fixed_sh.cc common_gamete.hpp
diploid_fixed_sh_ind_SOURCES=diploid_fixed_sh_ind.cc common_ind.hpp
diploid_fixed_sh_ind_lambda_SOURCES=diploid_fixed_sh_ind_lambda.cc common_ind.hpp
migsel_ind_SOURCES=migsel_ind.cc common_gamete.hpp
bneck_selection_ind_SOURCES=bneck_selection_ind.cc common_ind.hpp
diploid_ind_2locus_SOURCES=diploid_ind_2locus.cc common_ind.hpp
diploid_twopop_mig_SOURCES=diploid_twopop_mig.cc common_gamete.hpp
pfix_SOURCES=pfix.cc
RHH_SOURCES=RHH.cc
TFL2013_SOURCES=TFL2013.cc

AM_CPPFLAGS=-Wall -W -I.

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
LDADD=

