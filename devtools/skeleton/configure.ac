AC_PREREQ(2.59)

AC_INIT([FWDPPPACKAGE],[0.1.0],[FWDPPPROJECTURL])
AC_CONFIG_SRCDIR([src/FWDPPPACKAGE.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([enable])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile])

AC_ARG_ENABLE([standard],AS_HELP_STRING([--enable-standard],[Disable use of boost containers.]),[AC_DEFINE([USE_STANDARD_CONTAINERS],[1],[Prefer containers from namespace std when compiling])])

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

AC_CHECK_HEADER(boost/container/vector.hpp, [AC_DEFINE([HAVE_BOOST_VECTOR],[1],[Define if boost/container/vector.hpp exists])], [echo "boost/container/vector.hpp not found.  Will use std::vector for examples."])
AC_CHECK_HEADER(boost/container/list.hpp, [AC_DEFINE([HAVE_BOOST_LIST],[1],[Define if boost/container/list.hpp exists])], [echo "boost/container/list.hpp not found.  Will use std::list for examples."])
AC_CHECK_HEADER(boost/pool/pool_alloc.hpp, [AC_DEFINE([HAVE_BOOST_POOL_ALLOC],[1],[Define if boost/pool/pool_alloc.hpp exists])])
AC_CHECK_HEADER(boost/unordered_set.hpp,[AC_DEFINE([HAVE_BOOST_UNORDERED_SET],[1],[Define if boost/unordered_set.hpp exists])]) 
AC_CHECK_HEADER(boost/functional/hash.hpp, [AC_DEFINE([HAVE_BOOST_HASH],[1],[Define if boost/functional/hash.hpp exists])])
AC_CHECK_HEADER(Sequence/SimDataIO.hpp, , [AC_MSG_ERROR([Sequence/SimDataIO.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])

dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR(["zlib run time library not found";exit 1])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR(["gslcblas run-time library not found"])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR(["gsl run-time library not found"])])
dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[AC_MSG_ERROR([echo "libsequence runtime library does not seem to be here";exit 1])])
	  AC_CHECK_LIB(boost_system,main,,)
AC_LANG_RESTORE

AC_OUTPUT
