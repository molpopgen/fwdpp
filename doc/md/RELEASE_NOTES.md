# FWDPP RELEASE NOTES

## 0.3.7

* An declaration of KTfwd::recombine_gametes was fixed to match the implementation
* Minor re-organization of recombination code
* The details of KTfwd::ms_sample_separate were moved into KTfwd::internal.  This move allows access to functions that sample a specific set of diploids, which is handy in a lot of cases.

## 0.3.6

* The internal lookup table has been further refined, resulting in much faster simulations with selection.
* KTfwd::additive_diploid and KTfwd::multiplicative_diploid now have checks to prevent returning fitnesses < 0.

## 0.3.5

* Faster simulations with selection due to new data structure.  The gamete lookup tables introduced in 0.3.3 have been reworked.  This is not an API change, but rather an internal change that will be transparent to programmers using the library.
* "make check" now executes all unit tests, if and only if the boost unit testing library is present.
* Source code reorganization -- Issue #21
* Sugar code can now be included all at once.  This resolves Issue #23.  The convenience header is:

~~~{cpp}
#include <fwdpp/sugar.hpp>
~~~

* configure script now allows linking to the jemalloc library.  Preliminary benchmarks suggest that linking to tcmalloc gives better performance, but users may test for themselves:

~~~~{sh}
./configure --enable-jemalloc=yes
~~~~

* Travis-CI integration was added for the Git repo.
* Getting the Travis-CI integration to work resulted in changes to the unit tests and the experimental examples to improve their handling of boost-related things.
* unit/policyTests.cc was modified.  Previously, it was failing on GCC for the "wrong" reason.  The moves were actually happening, but the return value was being set at the wrong time.
* Various code cleanup

## 0.3.4

* LICENSE/COPYING files updated
* Tutorial improved for multi-locus simulations
* Issue #20 addressed throughout the library + examples + unit tests

##0.3.3

This release of __fwdpp__ includes major performance improvements.  The short version of the story is:

* Linked lists are traversed less often each generation compared to previous versions.
* An (extremely) expensive linear-time search of a linked list that previously occurred after _every_ recombination event has been eliminated, and replaced with a very fast log-time search of a lookup table.
* The containers required for storing the intermediate steps of recombination are now allocated at most once per replicate, and their capacity is adjusted as needed during the simulation.  Previously, these containers were allocated for every crossover event.
* The way that mutations are copied during recombination has been changed from a call to std::copy to the vector's insert member function, resulting in fewer reallocations during copying.

The result is much better scaling with large population size and recombination rates.

Major changes:

* Crossing over has been streamlined.  Unfortunately, this changes the library API.  However, run times improve substantially.
* The library internals now use different insertion methods during recombination.  We swtiched from copy(beg,end,back_inserter(x)) to x.insert(x.end(),beg,end), which results in less memory usage, and some run-time improvement for large simulations.
* API change: the simplification of metapopulation containers in 0.3.2 means that we can make the recombination policies required for such simulations the same as for single-deme simulations.
* The build system has been streamlined, with fewer dependencies and a greater emphasis on standard containers + Google's perftools.
* The various read/write functions for serializations have been renamed--I simply couldn't get the template declarations to auto-deduce the types.  I gave up on this one, and so there's an API change.  On the plus side, IO.hpp/IO.tcc are more readable now.

Minor changes:

* KTfwd::serialize now has a default constructor and a move constructor defined, allowing it to be a member of another class, which helps things in [foRward](http://github.com/molpopgen/foRward).
* KTfwd::serialize can now be used to serialize multiple records into a single buffer.  
* A "validation suite" has been added.  This is mostly for the developer to have an automatic way to check for problems.
* The package now installs a single binary called fwdppConfig.  This program is intended to be used in configure scripts to make checking for fwdpp's existence and/or version easier.  It takes a single option:

~~~{sh}
#Get the version of fwdpp installed on your system (or at least the one most readily visible in your user's environment).
fwdppConfig --version
~~~

## 0.3.2

This release make some tweaks that improve performance:

* All occurrences of boost::pool_allocator have been replaced by boost::fast_pool_allocator in the declarations of population types in the sugar layer (see @ref md_md_sugar).  This can result in >= 10% reductions in run times are recombination rates increase.
* The mechanics of adding mutations and recombining gametes have been streamlined slightly.
* A handful of cases where objects were copied instead of passed-by-reference were fixed.
* Issue #17 fixed
* Issue #14 fixed
* Issue #16 fixed

## 0.3.1

This release has many significant changes.  With the exception of the removal of KTfwd::tags::gamete_dependent, which was released in 0.3.0, and the change in the implementation of KTfwd::infsites::operator(), the library remains source-compatible with existing simulations.

* Multilocus, single-population version of KTfwd::sample_diploid now allows for inbreeding coefficient of 0 by default.
* Warnings from KTfwd::infsites suppressed w/compiler trick
* Extraneous gsl_rng free removed from diploid_ind.cc
* Issue #12 in examples/common_ind.hpp fixed
* Issue #13 appears to be fixed, too.  This was tricky, as it was not reproducible on my systems, but is now fixed on the system of the person who reported it.
* The library no longer requires that std::pair< T,T > are the types used for a diploid genotype (where the "T" are iterator types derived from gametes lists).  Now, custom diploid genotypes may be used, as long as the follow certain constraints (see @ref md_md_customdip for documentation).  Thanks to Jeremy Van Cleve for motivating this change, which should allow models involving space and social interactions, once some experimental features are worked out and become part of the main library (see below).
* extensions/ directory renamed python_examples
* Makefiles added for the boost.python examples.  These are not generated by ./configure, as compiling these is simply different from compiling the examples/unit tests
* Fitness policies for multilocus simulations now take iterators pointing to multilocus genotypes, in lieu of const-references to those genotypes.  This is done for the sake of API consistency across the library, but it does break source-level compatibility with existing simulations (although fixes will be easy).  Examples, unit tests, tutorials, etc., have been updated to account for this change.
* Mutation distpach methods implemented in 0.3.0 via KTfwd::tags::gamete_dependent are now dispatched using std::true_type and std::false_type
* Namespace KTfwd::experimental has been added as a safe place for experimenting with future library features.  The source files are in fwdpp/experimental. These headers do get installed with the library, but documentation will be minimal.
* Example programs based on experimental features are in experimental_examples.  These examples are "bleeding edge", and are not guaranteed to compile on all systems, especially those without complete working boost installations (including the compiled run-time libraries).
* experimental_examples/sex_limited_ind.cc is an example of using experimental library features to simulate separate sexes, sex-specific fitness effects, etc.
* Sugar typedefs in fwdpp/sugar/singlepop.hpp, fwdpp/sugar/metapop.hpp, and in fwdpp/sugar/multiloc.hpp are now compatible with custom diploid types.  The template alias now defaults to a pair of iterators pointing to gametes, which you may over-ride with your own type.
* The devtools setup script had some bugs fixed. (@ref md_md_devtools)
* KTfwd::sample_diploid has been modified to make fewer calls to random number generators for the common use case of no selfing.  This change results in N*generations fewer RNG calls, and can speed simulations by as much as ten percent. A side-effect is that results will differ from simulations compiled against earlier versions of the library.  A preprocessor macro/symbol (FWDPP_COMPAT_0_3_0) has been introduced.  Compilling programs with -DFWDPP_COMPAT_0_3_0 will result in the older/slower algorithm being used.
* A new document describing preprocessor symbols and their effect on fwdpp has been added (@ref md_md_preproc)
* All occurrences of NULL in the library have been replaced with the C++11 [nullptr](http://www.cplusplus.com/reference/cstddef/nullptr_t/).  The latter is safer in templates because the former can evaluate to 0 in some cases.
* All gamete-based code has been removed from the library.  If you need it for any reason, you're stuck with with 0.3.0 and earlier.
* Mutation models may now have various function signatures.  See @ref md_md_policies and the example programs for details.  The file fwdpp/tags/mutation_tags.hpp has been removed.  This was introduced in 0.3.0, and defined KTfwd::tags::gamete_dependent, which has been replaced with a newer, much more flexible method of policy dispatching.  See fwdpp/internal/mutation_internal.hpp for the details of how this was done.
* fwdpp/type_traits.hpp was added to the library, proving namespace KTfwd::traits.
* KTfwd::infsites no longer requires that a pointer to a mutation list be passed to the call operator
* KTfwd::mutation_writer and KTfwd::mutation_reader now support gzFiles.
* KTfwd::gzserialize and KTfwd::deserialize have been added to fwdpp/sugar/serialization.hpp

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
