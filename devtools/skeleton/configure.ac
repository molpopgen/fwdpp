AC_PREREQ(2.59)

dnl Default package version is 0.1.0, which is what KRT typically
dnl uses to correspond to an "initial published version"
dnl accompanying a paper
AC_INIT([FWDPPPACKAGE],[0.1.0],[FWDPPPROJECTURL])
AC_CONFIG_SRCDIR([src/FWDPPPACKAGE.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

FWDPPVERSION=`fwdppConfig --version`

dnl Check that fwdpp version is sufficient
dnl fwdpp versions >= 0.3.3 install a binary called fwdppConfig
dnl This program can be used to enable testing for specific versions of 
dnl the library.
dnl This if/else block is a simple method for doing that check.
dnl It will fail if:
dnl 1. fwdppConfig is not found in a user's $PATH
dnl 2. The return value of fwdppConfig --version fails the lexical string comparison done by "test"
if test "$FWDPPVERSION" \> "0.3.2"
then
	echo "fwdpp version $FWDPPVERSION detected"
else
	AC_MSG_ERROR([fwdpp >= 0.3.3 required, please install from http://github.com/molpopgen/fwdpp])
fi

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([disable])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile src/Makefile])

dnl enable use of boost list/vector/etc., if desired
AC_ARG_ENABLE([boost],AS_HELP_STRING([--enable-boost],[Enable use of boost containers when compiling examples]),[AC_DEFINE([USE_BOOST_CONTAINERS],[1],[Use containers from namespace boost when compiling])])

dnl tcmalloc from Google's perftools
AC_ARG_ENABLE([tcmalloc],AS_HELP_STRING([--enable-tcmalloc],[Enable linking to Google's tcmalloc library, if present.]))

dnl Do not compile with -DNDEBUG, so that (C-style) assertions are enables via <cassert>
AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

AC_CHECK_HEADER(fwdpp/diploid.hh,,[AC_MSG_ERROR([fwdpp/diploid.hh not found -- please install fwdpp (https://github.com/molpopgen/fwdpp)])])

dnl check for optional dependencies
if test x$enable_boost = xyes
then
AC_CHECK_HEADER(boost/container/vector.hpp, [AC_DEFINE([HAVE_BOOST_VECTOR],[1],[Define if boost/container/vector.hpp exists])], [echo "boost/container/vector.hpp not found.  Will use std::vector for examples."])
AC_CHECK_HEADER(boost/container/list.hpp, [AC_DEFINE([HAVE_BOOST_LIST],[1],[Define if boost/container/list.hpp exists])], [echo "boost/container/list.hpp not found.  Will use std::list for examples."])
AC_CHECK_HEADER(boost/pool/pool_alloc.hpp, [AC_DEFINE([HAVE_BOOST_POOL_ALLOC],[1],[Define if boost/pool/pool_alloc.hpp exists])])
AC_CHECK_HEADER(boost/unordered_set.hpp,[AC_DEFINE([HAVE_BOOST_UNORDERED_SET],[1],[Define if boost/unordered_set.hpp exists])]) 
AC_CHECK_HEADER(boost/functional/hash.hpp, [AC_DEFINE([HAVE_BOOST_HASH],[1],[Define if boost/functional/hash.hpp exists])])
fi

dnl required dependencies
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
	  dnl if we want boost containers, then we need to -lboost_system at link time
	  dnl Note that other useful boost libraries like interprocess, etc., also depend on the boost_system runtime lib,
	  dnl and so you may want to undo the conditional around this lib check.
	  if test x$enable_boost = xyes
	  then
		AC_CHECK_LIB(boost_system,main,,)
	  fi
AC_LANG_RESTORE

AC_OUTPUT
