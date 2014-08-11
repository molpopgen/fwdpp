AC_PREREQ(2.59)

AC_INIT([fwdpp], [0.2.5], [https://github.com/molpopgen/fwdpp])
AC_CONFIG_SRCDIR([fwdpp/diploid.hh])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config-h.in])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_C_CONST
AC_PROG_CXX

dnl enable c++11 support
AC_ARG_ENABLE(cpp11,
[   --enable-cpp11    Turn on debugging],
[case "${enableval}" in
   yes) AX_CXX_COMPILE_STDCXX_11([noext],[mandatory]) ;;
   no)  ;;
   *) AC_MSG_ERROR(bad value ${enableval} for --enable-cpp11) ;;
esac],[])


dnl AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

AM_MAINTAINER_MODE([enable])

dnl AC_PROG_LIBTOOL
AC_LANG(C++)
AC_CONFIG_FILES([Makefile fwdpp/Makefile examples/Makefile])

dnl AC_ARG_ENABLE(debug,
dnl [  --enable-debug    Turn on debugging],
dnl [case "${enableval}" in
dnl   yes) debug=true ;;
dnl   no)  debug=false ;;
dnl   *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
dnl esac],[debug=false])
dnl AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

dnl AC_ARG_ENABLE(profiling,
dnl [  --enable-profiling    Turn on profiling],
dnl [case "${enableval}" in
dnl   yes) profiling=true ;;
dnl   no)  profiling=false ;;
dnl   *) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
dnl esac],[profiling=false])
dnl AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)

dnl check for things that we need

AC_CHECK_HEADER(boost/bind.hpp,HAVE_BOOST=1 , [AC_MSG_ERROR([boost/bind.hpp not found.  boost libraries required.  Install boost (www.boost.org) ])])
AC_CHECK_HEADER(boost/container/vector.hpp, , [AC_MSG_ERROR([boost/container/vector.hpp not found.  Please install boost or upgrade your installation to a newer version])])
AC_CHECK_HEADER(boost/container/list.hpp, , [AC_MSG_ERROR([boost/container/list.hpp not found.  Please install boost or upgrade your installation to a newer version])])
AC_CHECK_HEADER(Sequence/SimDataIO.hpp,HAVE_LIBSEQ=1 , [AC_MSG_ERROR([Sequence/SimDataIO.hpp not found. libsequence (www.molpopgen.org) >= 1.7.8 required.])])
AC_CHECK_HEADER(gsl/gsl_rng.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(gsl/gsl_randist.h,HAVE_GSL=1,[AC_MSG_ERROR([gsl/gsl_rng.h not found.  Please install the GNU Scientific Library])])
AC_CHECK_HEADER(zlib.h,,[AC_MSG_ERROR([zlib.h not found.  zlib >= 1.2.5])])

dnl check for C run-time libraries
AC_CHECK_LIB([z],gzungetc,,[echo "zlib run time library not found";exit 1])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

dnl check for C++ run-time libraries
AC_LANG_SAVE
          AC_LANG_CPLUSPLUS
          AC_REQUIRE_CPP
          AC_CHECK_LIB(sequence, main,,[echo "libsequence runtime library does not seem to be here";exit 1])
AC_LANG_RESTORE

AC_OUTPUT