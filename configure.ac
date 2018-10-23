AC_PREREQ(2.59)

m4_define([FWDPP_MAJOR],0)
m4_define([FWDPP_MINOR],7)
m4_define([FWDPP_REVISION],0)
AC_INIT([fwdpp], [FWDPP_MAJOR.FWDPP_MINOR.FWDPP_REVISION], [https://github.com/molpopgen/fwdpp])
AC_DEFINE(FWDPP_MAJOR_VERSION,FWDPP_MAJOR,[Define major version number.])
AC_DEFINE(FWDPP_MINOR_VERSION,FWDPP_MINOR,[Define major version number.])
AC_DEFINE(FWDPP_REVISION_VERSION,FWDPP_REVISION,[Define major version number.])
AC_DEFINE(FWDPP_VERSION_STRING,"FWDPP_MAJOR.FWDPP_MINOR.FWDPP_REVISION",[Define the version string.])
AC_SUBST(FWDPP_MAJOR_VERSION,[FWDPP_MAJOR])
AC_SUBST(FWDPP_MINOR_VERSION,[FWDPP_MINOR])
AC_SUBST(FWDPP_REVISION_VERSION,[FWDPP_REVISION])
AC_SUBST(FWDPP_VERSION_STRING,[FWDPP_MAJOR"."FWDPP_MINOR"."FWDPP_REVISION])
AC_CONFIG_SRCDIR([fwdpp/diploid.hh])
AC_CONFIG_SRCDIR([examples/diploid_ind.cc])
AC_CONFIG_SRCDIR([testsuite/unit/fwdpp_unit_tests.cc])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

AM_MAINTAINER_MODE([disable])

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AC_LANG(C++)
AC_CONFIG_FILES([Makefile fwdpp/version.hpp fwdpp/Makefile fwdpp/io/Makefile fwdpp/algorithm/Makefile fwdpp/io/detail/Makefile
                 fwdpp/internal/Makefile fwdpp/tags/Makefile fwdpp/sugar/Makefile fwdpp/sugar/poptypes/Makefile
                 fwdpp/sugar/gsl/Makefile fwdpp/extensions/Makefile fwdpp/meta/Makefile
                 fwdpp/ts/Makefile fwdpp/ts/detail/Makefile examples/Makefile testsuite/Makefile src/Makefile doc/fwdpp.doxygen]) 

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)
dnl do we want tcmalloc?
AC_ARG_ENABLE([tcmalloc],AS_HELP_STRING([--enable-tcmalloc],[Enable linking to Google's tcmalloc library, if present.]))
dnl do we want jemalloc?
AC_ARG_ENABLE([jemalloc],AS_HELP_STRING([--enable-jemalloc],[Enable linking to the jemalloc library, if present.]))

dnl check for things that are required to compile all examples and/or use library at all
AC_CHECK_HEADER(Sequence/SimData.hpp, SIMDATAFOUND=1 , [AC_MSG_WARN([Sequence/SimData.hpp not found. Example programs will not be compiled.])])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
dnl AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5 is requred])])

AM_CONDITIONAL([HAVE_SIMDATA_HPP],test x$SIMDATAFOUND = x1)

dnl disabled in 0.4.4
dnl check for optional depdendencies
dnl if test x$enable_boost = xyes
dnl then
dnl AC_CHECK_HEADER(boost/container/vector.hpp,[AC_DEFINE([HAVE_BOOST_VECTOR],[1],[Is boost vector header found?])],[AC_MSG_ERROR([boost containers requested but boost/container/vector.hpp not foun])])
dnl AC_CHECK_HEADER(boost/container/list.hpp,[AC_DEFINE([HAVE_BOOST_LIST],[1],[Is boost list header found?])],[AC_MSG_ERROR([boost containers requested but boost/container/list.hpp not foun])])
dnl AC_CHECK_HEADER(boost/pool/pool_alloc.hpp,[AC_DEFINE([HAVE_BOOST_POOL_ALLOC],[1],[Is boost pool alloc header found?])],[AC_MSG_ERROR([boost containers requested but boost/container/bool_alloc.hpp not foun])])
dnl AC_CHECK_HEADER(boost/unordered_set.hpp,[AC_DEFINE([HAVE_BOOST_UNORDERED_SET],[1],[Is boost unordered set header found>])],[AC_MSG_ERROR([boost containers requested but boost/unordered_set.hpp not foun])])
dnl AC_CHECK_HEADER(boost/functional/hash.hpp,[AC_DEFINE([HAVE_BOOST_HASH],[1],[Is boost hash header found?])],[AC_MSG_ERROR([boost containers requested but boost/functional/hash.hpp not foun])])
dnl fi

AC_CHECK_HEADER(boost/test/unit_test.hpp, BUNITTEST=1,[echo "boost/test/unit_test.hpp not found. Unit tests will not be compiled."])

AM_CONDITIONAL([BUNIT_TEST_PRESENT], test x$BUNITTEST = x1)

AC_CHECK_HEADER(boost/program_options.hpp, BPROGOPT=1,[echo "boost/program_options.hpp not found. Some example programs will not be compiled."])

AM_CONDITIONAL([BPO_PRESENT], test x$BPROGOPT = x1)

dnl check for C run-time libraries
dnl AC_CHECK_LIB([z],gzungetc,,[AC_MSG_ERROR([zlib run time library not found])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])
dnl check for tcmalloc if desired
if test "x$enable_tcmalloc" = xyes
then
AC_CHECK_LIB([tcmalloc],[tc_cfree],,[AC_MSG_WARN([tcmalloc not found])])
fi
dnl check for jemalloc if desired
if test "x$enable_jemalloc" = xyes
then
AC_CHECK_LIB([jemalloc],[je_free],,[AC_MSG_WARN([jemalloc not found])])
fi
dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,HAVELIBSEQ=1,[AC_MSG_WARN([libsequence runtime library does not seem to be here.  Example programs will not be compiled])])
dnl	  AC_CHECK_LIB(boost_unit_test_framework,main,BUNITTESTLIB=1,
dnl                   [AC_MSG_WARN([boost unit testing framework not found. Test suite will not be compiled.])])
	  AM_CONDITIONAL([HAVE_LIBSEQ_RUNTIME], test x$HAVELIBSEQ = x1)
      AC_CHECK_LIB(boost_program_options,main,HAVEBPOLIB=1,[AC_MSG_WARN([boost_program_options runtime library does not
                                                            seem to be here.  Some example programs will not be
                                                            compiled])])
	  AM_CONDITIONAL([BPO_TEST_LIB_PRESENT], test x$HAVEBPOLIB = x1)
AC_LANG_RESTORE

dnl AC_CONFIG_HEADER(config.h)

AC_OUTPUT
