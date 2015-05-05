# FWDPP RELEASE NOTES

## 0.3.1

* Multilocus, single-population version of KTfwd::sample_diploid now allows for inbreeding coefficient of 0 by default.
* Warnings from KTfwd::infinite_sites suppressed w/compiler trick
* Extraneous gsl_rng free removed from diploid_ind.cc
* Issue #12 in examples/common_ind.hpp fixed
* The library no longer requires that std::pair< T,T > are the types used for a diploid genotype (where the "T" are iterator types derived from gametes lists).  Now, custom diploid genotypes may be used, as long as the follow certain constraints.  Thanks to Jeremy Van Cleve for motivating this change, which should allow models involving space and social interactions.
* extensions/ directory renamed python_examples
* Fitness policies for multilocus simulations now take iterators pointing to multilocus genotypes, in lieu of const-references to those genotypes.  This is done for the sake of API consistency across the library, but it does break source-level compatibility with existing simulations (although fixes will be easy).  Examples, unit tests, tutorials, etc., have been updated to account for this change.

## 0.3.0

* Allow for mutation policies that depdend on the current state of a gamete via KTfwd::tags::gamete_dependent.  Example HOC_ind.cc illustrates its use.
* Added dispatch tag to KTfwd::gamete_base, via KTfwd::tags::gamete_type_tag.  Currently, this feature is not used in the library.  Its addition affects nothing other than what happens during compilation, and it may be used in the future to add new features to the library, or it may be remvoed.
* Namespace KTfwd::sugar and \ref sugar functions added.  
* All _individual-based_ examples have been rewritten using the sugar code
* New tutorial added: \ref md_md_sugar
* New tutorial added: \ref md_md_datatypes
* The policies tutorial (\ref md_md_policies) has been rewritten
* Unit tests for sugar layer added

## 0.2.9

* Added new unit test (siteDepFitness) to check for the effect of issue #8.  Fortunately, fitness calculations were correct even with the bug, as the rest of the function did the right thing.
* Fixed issue #8, so that the code block in question will actually do what is intended.

## 0.2.8

* Fixed issue #5, which was a bug in migrating gametes in multi-population, individual-based simulations.
* Fixed issue #6, which was a bug in how parents were copied in multi-population, individual-based simulations.
* Fixed issue #7, which was a bug in when gametes were updated post-sampling in multi-population, individual-based simulations.
* Two new example programs added: migsel_split_ind and migsel_split_ind_list
* fwdpp/IO.tcc streamlined using C++11 "auto" instead of nasty typedefs

## 0.2.7 

* Versions ms_sample and ms_sample_separate that took containers of gametes as arguments are now compatible with individual-based simulations.
* The "devtools" stuff is greatly improved, and a new tutorial added on how to use it: @ref md_md_devtools

## 0.2.6 

* const mutation lists now passed to binary output routines
* The library internals now fully support C++11 move semantics
* Default policies now support C++11 "perfect forwarding"
* The function KTfwd::recombine_gametes, which is provided for recombination in individual-based simulations, has an overloaded version.  This new version takes a fixed set of positions representing recombination breakpoints, allowing more modeling flexibility and making unit testing easier.
* The mechanics of recombination in multilocus simulations moved to function multilocus_rec in namesapce KTfwd::fwdpp_internal.
* A unit test for the function KTfwd::fwdpp_internal::multilocus_rec was added.
* Added tutorial for multilocus simulation implmentation
* Minor cleanups to the build system
* Reorganization of fwdpp/IO.hpp and fwdpp/IO.tcc to reduce code duplication and provide output routines for multilocus simulations.  There are no longer separate read functions for gzFiles.  The necessary operations are handled automatically by overloads of template functions in namespace KTfwd::fwdpp_internal.
* Reorganization of fwdpp/sampling_functions.tcc: overloads of ms_sample are now implemented via calls to ms_sample_separate.  The results from the latter function are then merged using C++11 move semantics.  This reduces code redundancy, reduces the possible locations of bugs, and should keep efficiency about the same.

## 0.2.5

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

## 0.2.4

The following changes:

1. fwdpp/IO.hpp and fwdpp/IO.tcc were updated to include the ability to read/write binary-format data to/from gzip files.  The gzip compression uses [zlib](http://zlib.net).  Note that the combination of POSIX file locking + gzip output requires more work than "plain" binary output + file locking.  Please see the new examples.  Differences in how gztell vs. ftell mean that the index files generated are different.  See my [tutorial](https://github.com/molpopgen/BigDataFormats) on on "big data" file formats for more detail.
2.  examples/pfix.cc was added.  This estimates the fixation probability of a mutation subject to genic selection in a constant-sized population of N diploids.
3.  examples/diploid_gzbinaryIO.cc and examples/diploid_gzbinaryIO_ind.cc were added.  They illustrate the new functions mentioned in point 1.
4. The folder "test" was added.  This contains several Rmd ([R Markdown](http://rmarkdown.rstudio.com/)) documents that you can process either  in [R studio](http://www.rstudio.com/) or plain old [R](http://www.r-project.org).  In either case, you'll need [knitr](http://cran.r-project.org/web/packages/knitr/index.html) installed.
