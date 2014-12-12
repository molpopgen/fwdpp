#Tools to assist writing projects using fwdpp

Usage:

```{sh}
#Warning, arguments are not checkt
sh setup.sh directory projectname
cd directory/projectname
#edit configure.ac and replace FWDPPPACKAGE, FWDPPPACKAGE.cc, and FWDPPPROJECTURL with appropriate values, and 
#do the same for src/Makefile.am and rename the .cc file in src to match.
autoreconf -fi
autoheader
automake --add-missing -c
```



