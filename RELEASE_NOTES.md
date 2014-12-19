#0.2.5

The following changes:

1. fwdpp/diploid_functions_ind_multilocus.tcc was added.  This contains a more natural method of simulating mutiple partially-linked regions.  The header fwdpp/diploid_functions.hpp contains the relevant prototypes.  These routines are still in "beta" stage, and full support for simulations of partially-linked regions will be put off until a future version, likely 0.2.6 or later.
2. Header files have been reorganized.  fwdpp/diploid.hh is still the header to use.  The reorg has been for the developer's sanity.
3. C++11 support is now required.  A side-effect is that users may use what is currently the dev branch of libsequence, which requires c++11 and no longer requires boost.  However, simulations using fwdpp will still be a _lot_ faster with boost installed!
4. Internally, functions from namespace boost have been replaced with the namespace std equivalents provided by the C++11 standard.  
5. The internal library functions have been audited for performance.  This audit resulted in changing where some book-keeping functions were called, replacing linear-time search algorithms with logarithmic-time searches, and removing an unnecessary call to std::sort following crossover events.  The result of this audit was a big speedup of simulation run-times, but no change in simulation output (in other words, the same random number seed gives the same result for the same parameters for the same program in 0.2.4 and 0.2.5).
6.  autoconf stuff rewritten (configure.ac, Makefile.am, etc.)
7. Internally, the library has moved to a more "functional programming" style, emphasizing lambda expressions over function objects.  This change means that the library no longer needs to define function objects that are only used by the library "internals".  Such function objects have been removed, replaced with lambda expressions, and will no longer clutter the library documentation.
8. The code for recombination has been streamlined quite a bit.  Same algorithm (although noting the changes in point 5 above), but way fewer lines of code.
9. The devtools directory has been added.  It contains a script to setup packages using fwdpp.
10. Examples are now built via "make check" and not by default
11. Refactoring so that gamete- and individual- based methods share common code.
12. A new namespace (KTfwd::fwdpp_internal) now resides in fwdpp/internal.  This sub-namespace contains some of the more important inner-workings of the library that are used in several places.  Part of the goal of this sub-namespace is to separate the deterministic stuff from the stochastic stuff, in order for the latter to be unit-testable.
13. Boost's unit testing library is now used for unit testing.  To goal is to make sure that the stuff in fwdpp/internal all works.  The tests are compiled via "make check".
14. Use config.h to manage preprocessor stuff

#0.2.4

The following changes:

1. fwdpp/IO.hpp and fwdpp/IO.tcc were updated to include the ability to read/write binary-format data to/from gzip files.  The gzip compression uses [zlib](http://zlib.net).  Note that the combination of POSIX file locking + gzip output requires more work than "plain" binary output + file locking.  Please see the new examples.  Differences in how gztell vs. ftell mean that the index files generated are different.  See my [tutorial](https://github.com/molpopgen/BigDataFormats) on on "big data" file formats for more detail.
2.  examples/pfix.cc was added.  This estimates the fixation probability of a mutation subject to genic selection in a constant-sized population of N diploids.
3.  examples/diploid_gzbinaryIO.cc and examples/diploid_gzbinaryIO_ind.cc were added.  They illustrate the new functions mentioned in point 1.
4. The folder "test" was added.  This contains several Rmd ([R Markdown](http://rmarkdown.rstudio.com/)) documents that you can process either  in [R studio](http://www.rstudio.com/) or plain old [R](http://www.r-project.org).  In either case, you'll need [knitr](http://cran.r-project.org/web/packages/knitr/index.html) installed.
