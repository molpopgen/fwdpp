AUTOMAKE_OPTIONS = foreign
SUBDIRS=fwdpp examples testsuite src
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST=doc/fwdpp.doxygen examples/*.cc examples/*.hpp LICENSE 

clean-local:
	find . -name '*.gcov' | xargs rm -f
	rm -f testsuite/*.info
	rm -rf fwdpp_coverage
