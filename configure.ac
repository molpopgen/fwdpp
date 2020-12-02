AC_PREREQ(2.59)

m4_define([FWDPP_MAJOR],0)
m4_define([FWDPP_MINOR],10)
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

AC_ARG_ENABLE([cpp17], AS_HELP_STRING([--enable-cpp17], [Compile using -std=c++17]))
if test "x$enable_cpp17" = x
then
AX_CXX_COMPILE_STDCXX([14],[noext],[mandatory])
else
AX_CXX_COMPILE_STDCXX([17],[noext],[mandatory])
fi

AC_LANG(C++)
AC_CONFIG_FILES([Makefile fwdpp/version.hpp fwdpp/Makefile
                 fwdpp/fundamental_types/Makefile fwdpp/io/Makefile
                 fwdpp/algorithm/Makefile fwdpp/io/detail/Makefile
                 fwdpp/types/Makefile
                 fwdpp/internal/Makefile fwdpp/tags/Makefile
                 fwdpp/sugar/Makefile fwdpp/poptypes/Makefile
                 fwdpp/gsl/Makefile fwdpp/extensions/Makefile
                 fwdpp/meta/Makefile fwdpp/simfunctions/Makefile
                 fwdpp/ts/Makefile fwdpp/ts/detail/Makefile
                 fwdpp/ts/marginal_tree_functions/Makefile
                 fwdpp/ts/simplification/Makefile
                 fwdpp/ts/recording/Makefile
				 fwdpp/util/Makefile
				 fwdpp/genetic_map/Makefile
				 examples/Makefile testsuite/Makefile
				 src/Makefile
				 .circleci/config.yml
				 ]) 

AC_ARG_ENABLE([debug],AS_HELP_STRING([--enable-debug],[Compile examples with run-time assert macro enabled.]))
AM_CONDITIONAL([DEBUG],test x$enable_debug = xyes)

dnl check for things that are required to compile all examples and/or use library at all
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])

AM_CONDITIONAL([HAVE_SIMDATA_HPP],test x$SIMDATAFOUND = x1)

AC_CHECK_HEADER(boost/test/unit_test.hpp, BUNITTEST=1,[echo "boost/test/unit_test.hpp not found. Unit tests will not be compiled."])

AM_CONDITIONAL([BUNIT_TEST_PRESENT], test x$BUNITTEST = x1)

AC_CHECK_HEADER(boost/program_options.hpp, BPROGOPT=1,[echo "boost/program_options.hpp not found. Some example programs will not be compiled."])

AM_CONDITIONAL([BPO_PRESENT], test x$BPROGOPT = x1)

dnl check for C run-time libraries
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_ERROR([gslcblas run-time library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([gsl run-time library not found])])

dnl check for C++ run-time libraries
AC_LANG_SAVE
      AC_LANG_CPLUSPLUS
      AC_REQUIRE_CPP
      AC_CHECK_LIB(boost_program_options,main,HAVEBPOLIB=1,[AC_MSG_WARN([boost_program_options runtime library does not
                                                            seem to be here.  Some example programs will not be
                                                            compiled])])
	  AM_CONDITIONAL([BPO_TEST_LIB_PRESENT], test x$HAVEBPOLIB = x1)
      AC_CHECK_LIB(boost_unit_test_framework,main,HAVEBOOSTTESTRUNTIME=1,[AC_MSG_WARN([boost_unit_test_framework runtime library does not
                                                            seem to be here.  Unit tests will not be compiled])])
      AM_CONDITIONAL([BUNIT_TEST_RUNTIME_PRESENT], test x$BUNITTEST = x1)
AC_LANG_RESTORE

dnl AC_CONFIG_HEADER(config.h)

AC_OUTPUT
