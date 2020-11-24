AUTOMAKE_OPTIONS = foreign
SUBDIRS=fwdpp examples testsuite src
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST=examples/*.cc examples/*.hpp LICENSE \
		   doc/conf.py \
		   doc/fwdpp.bib \
		   doc/index.rst \
		   doc/requirements.txt \
		   doc/Makefile \
		   doc/rst/*.rst

clean-local:
	find . -name '*.gcov' | xargs rm -f
	rm -f testsuite/*.info
	rm -rf fwdpp_coverage
