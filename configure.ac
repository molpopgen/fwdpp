AC_PREREQ(2.59)

AC_INIT([fwdpp], [0.2.5], [https://github.com/molpopgen/fwdpp])
AC_CONFIG_SRCDIR([fwdpp/diploid.hh])
AC_CONFIG_SRCDIR([examples/diploid_ind.cc])
AC_CONFIG_SRCDIR([unit/mutateTest.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config-h.in])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([disable])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile fwdpp/Makefile fwdpp/internal/Makefile examples/Makefile unit/Makefile])

AC_ARG_ENABLE([standard],AS_HELP_STRING([--enable-standard],[Disable use of boost containers when compiling examples]))
AM_CONDITIONAL([USE_STANDARD],test x$enable_standard = xyes)

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

dnl check for things that we might need
AC_CHECK_HEADER(boost/container/vector.hpp, BVECTOR=1, [echo "boost/container/vector.hpp not found.  Will use std::vector for examples."])
AC_CHECK_HEADER(boost/container/list.hpp, BLIST=1, [echo "boost/container/list.hpp not found.  Will use std::list for examples."])
AC_CHECK_HEADER(boost/test/unit_test.hpp, BUNITTEST=1,[echo "boost/test/unit_test.hpp not found. Unit tests will not be compiled."])
AC_CHECK_HEADER(boost/pool/pool_alloc.hpp, BPOOLALLOC=1)
AC_CHECK_HEADER(boost/unordered_set.hpp, BUNORDEREDSET=1)
AC_CHECK_HEADER(boost/functional/hash.hpp, BHASH=1)
AC_CHECK_HEADER(Sequence/SimDataIO.hpp, , [AC_MSG_ERROR([Sequence/SimDataIO.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])
AM_CONDITIONAL([BVEC_PRESENT],test x$BVECTOR = x1)
AM_CONDITIONAL([BLIST_PRESENT],test x$BLIST = x1)
AM_CONDITIONAL([BUNIT_TEST_PRESENT], test x$BUNITTEST = x1)
AM_CONDITIONAL([BPOOL_ALLOC_PRESENT], test x$BPOOLALLOC = x1)
AM_CONDITIONAL([BUO_SET_PRESENT], test x$BUNORDEREDSET = x1)
AM_CONDITIONAL([BHASH_PRESENT], test x$BHASH = x1)

dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR(["zlib run time library not found";exit 1])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR(["gslcblas run-time library not found"])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR(["gsl run-time library not found"])])
dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([eecho "libsequence runtime library does not seem to be here";exit 1])])
	  AC_CHECK_LIB(boost_system,main,,)
AC_LANG_RESTORE

dnl AM_CONFIG_HEADER(config.h)

AC_OUTPUT
