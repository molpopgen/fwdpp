AC_PREREQ(2.59)

AC_INIT([FWDPPPACKAGE],[0.1.0],[FWDPPPROJECTURL])
AC_CONFIG_SRCDIR([src/qRHH_single.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config-h.in])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([enable])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile])

AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

dnl check for things that we might need
AC_CHECK_HEADER(Sequence/SimDataIO.hpp,HAVE_LIBSEQ=1 , [AC_MSG_ERROR([Sequence/SimDataIO.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_CHECK_HEADER(fwdpp/diploid.hh,HAVE_FWDPP=1 , [AC_MSG_ERROR([fwdpp/diploid.hh not found. fwdpp >= 0.2.5 required])])
AC_CHECK_HEADER(fwdpp/rec_gamete_updater.hpp,HAVE_FWDPP_OK=1 , [AC_MSG_ERROR([fwdpp/rec_gamete_updater.hpp not found. fwdpp >= 0.2.5 required])])
AC_CHECK_HEADER(boost/container/vector.hpp, BVECTOR=1, [AC_MSG_WARN([echo "boost/container/vector.hpp not found.  Will use std::vector to compile."])])
AC_CHECK_HEADER(boost/container/list.hpp, BLIST=1, [echo "boost/container/list.hpp not found.  Will use std::list to compile."])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])
AM_CONDITIONAL([BVEC_PRESENT],test x$BVECTOR$enable_boost = x1)
AM_CONDITIONAL([BLIST_PRESENT],test x$BLIST = x1)

dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR(["zlib run time library not found";exit 1])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR(["gslcblas run-time library not found"])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR(["gsl run-time library not found"])])

dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([eecho "libsequence runtime library does not seem to be here";exit 1])])
AC_LANG_RESTORE

dnl AM_CONFIG_HEADER(config.h)

AC_OUTPUT