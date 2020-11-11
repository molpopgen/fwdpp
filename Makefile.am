AUTOMAKE_OPTIONS = foreign
SUBDIRS=fwdpp examples testsuite src
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST=doc/fwdpp.doxygen examples/*.cc examples/*.hpp LICENSE 

coverage-local: check testsuite/coverage.sh
	sh testsuite/coverage.sh

clean-local:
	find . -name '*.gcov' | xargs rm -f
	rm -rf fwdpp_coverage
