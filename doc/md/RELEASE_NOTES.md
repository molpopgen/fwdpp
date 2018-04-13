# FWDPP RELEASE NOTES 
For a list of planned features, etc., see the issues list on [GitHub](https://github.com/molpopgen/fwdpp/issues).
Issues that are tagged 'performance' or 'enhancement' reflect future plans for the library. I will probably not put
milestones (target version numbers for these features to go live) because that is not realistic given how I work.

## 0.6.0

This is a major release, breaking API compatibility in many areas.  These changes were necessary to make the library
even more flexible, but also to undo several poor design choices that we've been living with for far too long.  This
release can be seen as a big push towards a 1.0 release.  

The result of these changes is a much smaller library that is safer and can do more.

The following bug fixes are included in this release.  None affected
[fwdpy11's](http://fwdpy11.readthedocs.io/en/latest/) official releases, either because that
package uses custom replacements of fwdpp functions and/or they were caught before a release was made or are in features
so new that no one has affected by them:

* A bug in which selected fixations were not recorded by fwdpp::update_mutations_n was fixed.
  [PR119](https://github.com/molpopgen/fwdpp/pull/119)
* A bug in repeatedly recording fixation data with removeFixed == false has been fixed. 
  [PR71](https://github.com/molpopgen/fwdpp/pull/71)
* API bugs in multilocus sampling routines have been fixed. 
  [PR70](https://github.com/molpopgen/fwdpp/pull/70)
* The "perfect-forwarding" constructor of fwdpp::popbase now is indeed perfectly-forwarding.[PR65](https://github.com/molpopgen/fwdpp/pull/65).
  This fixes a bug introduced in [PR56](https://github.com/molpopgen/fwdpp/pull/56).

In fwdpp 0.5.7, [PR56](https://github.com/molpopgen/fwdpp/pull/56) had a few issues resolved in this release:

* Do not require mutation keys to be sorted by position in extinct gametes. [PR66](https://github.com/molpopgen/fwdpp/pull/66)  
* Fixed an error in a "perfect forwarding" constructor introduced in 0.5.7. [PR65](https://github.com/molpopgen/fwdpp/pull/65) 
* When creating populations based on user data, 0.5.7 failed to enforce
  proper sorting of mutation keys in gametes.  This issue was fixed via [PR64](https://github.com/molpopgen/fwdpp/pull/64).
    
The following changes have been made to the library itself:

* @ref custom_mutation.cc, @ref custom_diploid.cc and @ref juvenile_migration.cc were added to examples. [PR124](https://github.com/molpopgen/fwdpp/pull/124).
* All use of std::bind is replaced with lambdas. [PR126](https://github.com/molpopgen/fwdpp/pull/126).
* Some mutation types were removed from fwdpp/sugar. [PR125](https://github.com/molpopgen/fwdpp/pull/125).
* The back end for specifying final values from fwdpp::additive_diploid and fwdpp::multiplicative_diploid have been
  streamlined. [PR121](https://github.com/molpopgen/fwdpp/pull/121), which changes behavior introduced in fwdpp 0.5.6
  in response to [Issue 49](https://github.com/molpopgen/fwdpp/issues/49)
* The namespace has been changed from KTfwd to fwdpp. [PR88](https://github.com/molpopgen/fwdpp/pull/88)
* Serialization code has been generalized to depend on template specializations. [PR90](https://github.com/molpopgen/fwdpp/pull/90) and [PR108](https://github.com/molpopgen/fwdpp/pull/108)
* The struct KTfwd::infsites was removed.  For the mutation type fwdpp::popgenmut, fwdpp::infsites_popgenmut was added
  as a mutation function.  [PR120](https://github.com/molpopgen/fwdpp/pull/120)
* The storage of a gsl_rng pointer in fwdpp::poisson_xover and fwdpp::general_rec_variation was refactored out. [PR118](https://github.com/molpopgen/fwdpp/pull/118)
* The types fwdpp::extensions::discrete_rec_model and fwdpp::extensions::discrete_mut_model were refactored to be much
  more general, and independent of a particular mutation type. [PR116](https://github.com/molpopgen/fwdpp/pull/116) and [PR113](https://github.com/molpopgen/fwdpp/pull/113)
* Metapopulation objects were removed, and single-/multi- locus population types are renamed fwdpp::slocuspop and  fwdpp::mlocuspop, respectively.  Code for demographic operations was also removed, as it is not needed. 
  [PR110](https://github.com/molpopgen/fwdpp/pull/110)
* The class fwdpp::popbase is no longer dependent on a particular ploidy. The diploid-ness was pushed into fwdpp::slocuspop and fwdpp::mlocuspop (as well as the more complex template classes fwdpp::sugar::slocuspop and fwdpp::sugar::mlocuspop). [PR114](https://github.com/molpopgen/fwdpp/pull/114)
* Namespace KTfwd::experimental was removed.  [PR86](https://github.com/molpopgen/fwdpp/pull/86)
* Diploid-dependent mutation models are now supported. [PR84](https://github.com/molpopgen/fwdpp/pull/84)
* Removed API features deprecated in 0.5.7.  [PR83](https://github.com/molpopgen/fwdpp/pull/83)
* Several parts of namespace fwdpp::traits were streamlined. [PR82](https://github.com/molpopgen/fwdpp/pull/82)
* The "scaling" parameter of fwdpp::additive_diploid and fwdpp::multiplicative_diploid are now class data instead of
  bound parameters.  [PR79](https://github.com/molpopgen/fwdpp/pull/79)
* fwdpp::general_rec_variation added. [PR77](https://github.com/molpopgen/fwdpp/pull/77)
* Recombination callbacks taking no arguments are supported. [PR75](https://github.com/molpopgen/fwdpp/pull/75)
* fwdpp::mlocuspop::locus_boundaries is no longer populated with default values, which was dangerous!! 
  [PR69](https://github.com/molpopgen/fwdpp/pull/69)
* fwdpp::generalmut and fwdpp::generalmut_vec now take less memory.
  [PR67](https://github.com/molpopgen/fwdpp/pull/67)
* Sampling routines in fwdpp/sugar/sampling.hpp now guard against the possibility of having fixations repeated in
  samples. For example, this change will affect simulations where selected fixations are retained because they affect trait values. [PR71](https://github.com/molpopgen/fwdpp/pull/71) 
* Issue [PR69](https://github.com/molpopgen/fwdpp/issues/69) [PR70](https://github.com/molpopgen/fwdpp/pull/70) 
* Fixed template errors when sampling from multilocus populations using functions in fwdpp/sugar/sampling.hpp.  These
  errors resulted in failure to compile and thus did not affect anyone's results.  [PR70](https://github.com/molpopgen/fwdpp/pull/70) 


## 0.5.7

* License change from GPL2 to GPL3
* The unit test suite is now compiled with warnings enabled.  As a result, several warnings were silenced, including one
  introduced in [PR55](https://github.com/molpopgen/fwdpp/pull/56) that would have mattered!.  [PR63](https://github.com/molpopgen/fwdpp/pull/63)
* API to KTfwd::sample_diploid was updated to removed an unused type as a result of [PR54](https://github.com/molpopgen/fwdpp/pull/54). [PR63](https://github.com/molpopgen/fwdpp/pull/63)
* Gametes can now be constructed from tuples. [PR62](https://github.com/molpopgen/fwdpp/pull/62)
* KTfwd::popgenmut, KTfwd::generalmut, and KTfwd::generalmut_vec may now be constructed from tuples. [PR59](https://github.com/molpopgen/fwdpp/pull/59)
* Travis build system now skips OS X builds, tries more GCC versions as well as C++11 and C++14.
  [PR60](https://github.com/molpopgen/fwdpp/pull/60) and [PR61](https://github.com/molpopgen/fwdpp/pull/61).
  We are skipping OS X because there have been stability problems with the Travis CI infrastructure for that platform.
  We would liked to have included Linux/clang tests for this release, but hit road blocks that we cannot reproduce on our
  local machines.
* Population objects from sugar layer can now be constructed with pre-calculated diploids, gametes, and mutations.
  [PR56](https://github.com/molpopgen/fwdpp/pull/56).
* Mutation and recombination are now merged into one path.  See fwdpp/mutate_recombine.hpp.  The entry points into the
  old API are marked deprecated.  This addresses issue [Issue 54](https://github.com/molpopgen/fwdpp/issues/54) via pull request [PR55](https://github.com/molpopgen/fwdpp/pull/55).
* Travis builds using miniconda now use -Wl,-rpath when linking to dependencies, solving run-time link errors during "make check" on OS X.

## 0.5.6

This release breaks API compatibility.

* Issue #51 resolved.
* API issue #50 resolved.  The KTfwd::GSLrng_t is no longer copy-constructible.
* API issue #48 is addressed.  This is an API change breaking backwards compatilibity for multi-locus/region simulations. The new API allows more flexibility in modeling interlocus/region crossovers.  The API for KTfwd::sample_diploid is changed for multi-locus/region sims.  Two functions, KTfwd::make_poisson_interlocus_rec and KTfwd::make_binomial_interlocus_rec, return vectors of callbacks bound to the two new structs. [commit](https://github.com/molpopgen/fwdpp/commit/8ee950e7f315434284164e50c0f09b1e52a4c40c)
* API issue #49 is addressed.  The changes maintain compile-time API compatibility with previous library versions.
* Refactored KTfwd::extensions::discrete_rec_model and KTfwd::extensions::discrete_mut_model to use PIMPL idiom and be
  default-constructible. [commit](https://github.com/molpopgen/fwdpp/commit/9edcb8ca0da2dac5d04f066fbc0f26e3b7776c16)
* Extend concept of regions to multi-locus sims via extensions::bind_vec__drm and extensions::bind_vec::dmm. [commit](https://github.com/molpopgen/fwdpp/commit/da1b47b661679c80530b3ed477107f9fadf33e25)
* extensions::discrete_mut_model::make_mut now takes a pointer to the generation, allowing a single point of binding rather than binding each generation. [commit](https://github.com/molpopgen/fwdpp/commit/da1b47b661679c80530b3ed477107f9fadf33e25)
* Exceptions in namespace KTfwd::extensions changed from std::runtime_error to std::invalid_argument where appropriate.
* Fix bug in examples/HOC_ind.cc affecting (improper) recycling of fixations during the simulation. [commit](https://github.com/molpopgen/fwdpp/commit/559e7db4cefe6c444584c4a51587bd315e35cbb9)
* KTfwd::data_matrix is now based on std:int8_t instead of char.

## 0.5.5

* Added KTfwd::sugar::multiloc::locus_boundaries and simplified API in fwdpp/sugar/sampling.hpp, resolving #47.
* Refactor namespace KTfwd::traits, resolving #46.
* Streamline API for fitness functions, resolving #45.
* Resolving #45 allows us to simplify how custom diploid types are defined.  See @ref md_md_customdip
* Test suite fixtures improved
* All source files reformatted using clang-format
* Build system generation of package version numbers improved
* fwdpp/version.hpp added and is auto-generated by build system
* fwdpp/sugar/matrix.hpp now included by fwdpp/sugar/sugar.hpp
* Fix forwarding of constructor arguments during mutation recycling.
* python_examples/fwdpp_pybind11.cc was updated to account for population type class hierarchy added in 0.5.3.

## 0.5.4

* fwdpp/sugar/matrix.hpp was added, providing functions for returning diploid haplotype/genotype data as a 1-d array.
* Functions in fwdpy/sugar/sampling.hpp are now more flexible with the integral types used to specify specific sets of diploids.
* Fixed API bugs in KTfwd::sample_separate for the case of multi-locus/region simulations.  These errors resulted in compilation failure, and therefore it was impossible for results to have been affected.

## 0.5.3 

* Streamline internal details of crossing over
* Fix issue #43
* KTfwd::change_neutral no longer tries to update extinct gametes.
* Single-argument constructors marked explicit for types inheriting from KTfwd::popbase
* Added an additional overload of KTfwd::infsites::operator()
* The integer type stored by a gamete and used to index mutations was changed from std::size_t to std::uint32_t.  This change halves RAM use and has no other side-effects other than limiting the number of possible mutations in a simulation to 2^32, which is too many to store on a typical cluster node anyways.

There are also several changes to the build setup and Travis CI:

* Travis CI is now based on miniconda, taking advantage of bioconda/libsequence.
* Travis CI now only tests GCC on OS X, as bioconda/libsequence is built with that compiler
* "make" now makes fwdppConfig, all examples, and all unit tests, provided that dependencies are present.  The examples and unit tests are not installed.
* "make check" now runs the unit tests

## 0.5.2

* Documentation updates, finally!  The tutorials, etc., have been brought up to date.
* Missing include of cassert added to fwdpp/internal/recycle.hpp
* Serialization code streamlined.  Lots of redundant code was removed.  The biggest changes are that std::runtime_error can be thrown from low-level functions.  Also, KTfwd::serialize can now work with any stream type whose public interface is compatible with std::istream, and KTfwd::gzserialize is simply a convenience wrapper (it still works via an in-memory serialization).
* Test suite refactored.  The code has moved from unit/ to testsuite/, and attempts to better separate unit tests from
  integration tests.  Further, fixtures are used to improve code reuse in testing.
* operator== added to KTfwd::generalmut and KTfwd::generalmut_vec

## 0.5.1

This is a small bigfix release.  There were performance improvements planned for this release,
but we will hold off in favor of getting this fix out.

* Fix for #41 in v0.5.0 did not get applied to the experimental version of sample_diploid for multi-locus sims.  That has been fixed.

## 0.5.0

* Example file examples/K_linked_regions_multilocus.cc added
* Streamlined KTfwd::fwdpp_internal::multilocus_rec_mut
* The sugar types KTfwd::singlepop, KTfwd::metapop, and KTfwd::multiloc were refactored to inerit from KTfwd::sugar::popbase.
* Fixed error in KTfwd::multiloc where gametes were initialized with the incorrect count
* KTfwd::gamete_data_sane_multiloc was added to fwdpp/debug.hpp
* Issue #41 fixed. This issue affected simulations using the multi-locus API, and all simulations using that API need to be rerun. Sorry.
* The experimental API to sample_diploid was made more flexible via a new header file, fwdpp/experimental/dispatch.hpp.  This addition allows better fine-tuning of "rules" classes
* The experimental API now takes rvalue refrence (&&) instead of const reference (const &) for rules classes.  This allows the rules to be written and be more idiomatic, avoiding use of mutable variables.
* Add overload of KTfwd::haplotype_dependent_fitness::operator() for custom diploids

## 0.4.9

* Doxygen file changed so that all library source is browsable.
* Single-deme version of sample_diploid is commented in detail.
* Added KTfwd::add_mutation, KTfwd::add_mutations and unit test  unit/test_sugar_add_mutation.cc.  These new functions allow the addition of mutations to diploids in a non-random way, which also means you can fill in a population from external data (e.g., something in a file).  These features resolve Issue #28, albeit at a low level.
* Member KTfwd::mutation_base::xtra squeezed into unused space in this type.  No extra RAM used, and programs may assign values to that type to represent "stuff", whatever that is.
* KTfwd::extensions::discrete_mut_model got a new constructor allowing KTfwd::extensions::discrete_mut_model::make_mut to assign values to KTfwd::mutation_base::xtra.
* Added KTfwd::change_neutral, which allows simulations to update the value of KTfwd::mutation_base::neutral and correctly update storage of the affected mutation in all gametes.
* #39 fixed
* #40 fixed

## 0.4.8

* #38 fixed
* Updates to unit tests
* Added and example of wrapping fwdpp using [pybind11](https://github.com/pybind/pybind11)

## 0.4.7

Thanks to Alexander Nater for pointing out issue #36, which lead to #37 being discovered, too.

* remove include<iostream> from fwdpp/internal/recombination_common.hpp
* Isssue #36 fixed
* Issue #37 fixed
* KTfwd::GSLrng_t has improved copy constructor and is now move-constructible.
* The library is now const-correct vis-a-vis gsl_rng *.  This is an API change.
* "debug mode" (compiling _without_ -DNDEBUG) makes more tests about samples from populations being sorted by position
* More extensive checking done by functions in fwdpp/debug.hpp
* Improvements (and bug fixes) to namespace KTfwd::traits.  Unit tests added for this namespace.
* Fixed templates for serializing mutation type KTfwd::generalmut.  Attempting to compile using icc revealed the error.
* KTfwd::fwdpp_internal::gamete_cleaner has a new, and generally much faster, implementation.

## 0.4.6

* Issue #34 fixed
* Issue #35 fixed.  Thanks to Alexander Nater for catching this.
* KTfwd::extensions::discrete_rec_model no longer generates empty lookup tables when there are no regions
* KTfwd::extensions::discrete_rec_model and KTfwd::discrete_mut_model may now throw exceptions from their constructors if input data are incorrect

## 0.4.5

* Issues #31 and #32 fixed
The file fwdpp/internal/gamete_lookup_table.hpp was removed and the library updated to stop using this method to check for gamete uniqueness.  The result is a 20% reduction in run times and a slight reduction in peak memory use.  This change was enabled by the changes introduced in 0.4.4, and comes with no change in output.
* Partial loop unrolling and branch removal from KTfwd::fwdpp_internal::recombine_gametes
* Faster fitness calculations via the removal of if statements from a for loop in KTfwd::site_dependent_fitness
* KTfwd::fwdpp_internal::add_new_mutation was shortened for clarity
* KTfwd::metapop objects now copy- and move- constructible from KTfwd::singlepop objects
* Support for demographic events via the low-level functions KTfwd::copy_deme, KTfwd::merge_demes, KTfwd::remove_deme, KTfwd::swap_demes, KTfwd::split_deme, and KTfwd::admix_demes.
* Higher-level support for demographic models via the sugar-layer functions KTfwd::copy_pop, KTfwd::merge_pop, KTfwd::remove_pop, KTfwd::swap_pops, KTfwd::split_pop, and KTfwd::admix_pops.
* New unit tests added and old ones refined

## 0.4.4

TL;DR Big performance improvements due to better handling of objects in memory. Yay!  API changes. Boo!

__The library documentation is likely to be out of date until a future release__

* Vastly improved management of object lifetimes:
  * Extinct mutations/gametes are no longer removed each generation.
  * FIFO queues are constructed each generation in order to "recycle" those objects as new mutations/gametes.
    * The queue is implemented as std::queue<mlist::iterator> or std::queue<glist::iterator> for mutations and gametes, respectively.
    * These queues result in big performance improvements for "big" simulations, at the cost of breaking API compatibility.  The performance improvement is both reduced run time and reduced memory usage.
    * In order to accomodate the "recycling", mutation models must now return iterators derived from the mutation list, rather than mutation types themselves.
    * Additionally, the internal mutation/recombination functions must take non-const references to "recycling bins", which are the FIFO queues.
  * Mutation type data members (objects inheriting from KTfwd::mutation_base) are no longer const.  This is required in order to enable the "recycling".
  * In order to take advantage of this feature, extinct mutations must not be removed each generation.  The function KTfwd::update_mutations should be removed instead of KTfwd::remove_lost, KTfwd::remove_fixed_lost, etc.  The latter functions are marked as deprecated.
  * The details of recycling are handled by KTfwd::fwdpp_internal::make_mut_queue, KTfwd::fwdpp_internal::make_gamete_queue, KTfwd::fwdpp_internal::recycle_gamete, and the very cool variadic template function KTfwd::fwdpp_internal::recycle_mutation_helper.  These types/functions are found in fwdpp/internal/recycling.hpp.

* Improved memory managament:
  * "Recycling" means that the linked lists used are much more constant in memory with respect to allocations.  Thus, it makes sense to allocate objects prior to simulation.
  * The functions KTfwd::add_elements in the main library and KTfwd::add_recyclable in the sugar sub-library allow the addition of recylable objects.  See example programs for how to use these functions.

* Other API changes:
  * As a result of the changes described above, the data structures of the library were changed from linked lists of iterators to vectors of std::size_t.  This is a radical change to the data layout, and results in further speedups.
  * KTfwd::sample_diploid has been streamlined in light of common defaults
  * Recombination policies have been simplified.  See KTfwd::poisson_xover for an example.

* Behavior changes:
  * As a result of object recycling, data structures (mutation and gamete lists, specifically) at the end of a simulation contain both extant and extinct objects.
* Implementation changes:
  * Serialized populations may contain extinct mutations and gametes.
  * A lot of redundant code has been replaced with function calls. This should help prevent bugs like Issue #27, which was due to a botched copy-paste, which happens when code is written in a hurry...
  * Simpler dispatch method for mutation models (KTfwd::fwdpp_internal::mutation_model_dispatcher)

* Improved "type traits" sub-library (namespace KTfwd::traits)
* Deprecated functions/objects removed from library
* Bug fixes:
  * Issue #30 fixed regarding serialization of KTfwd::generalmut and KTfwd::generalmut_vec
  * A bug in KTfwd::ms_sample and KTfwd::ms_sample separate was identified and fixed.  The bug only applied to multilocus simulations.  The bug was that the first 'n' chromosomes were sampled, rather than 'n' randomly-chosen chromosomes.  The effect of the bug is that the first sample is truly random (as chromosomes are not sorted in any meaningful way), but a second sample would overlap with the first, etc.

## 0.4.3

* Fix for issue #29
* Keyword 'mutable' replaced with 'const' throughout library
* KTfwd::extensions::gaussian now uses ziggurat method
* Types declared in fwdpp/extensions/callbacks.hpp now have const member data.  A unit test was added as a check on the API of this file.
* New mutation types added: KTfwd::generalmut and KTfwd::generalmut_vec, in fwdpp/sugar/generalmut.hpp.  A new unit test file goes along with it.

## 0.4.2

* Fixed error in definition of KTfwd::metapop_serialized and KTfwd::multiloc_serialized.

## 0.4.1

* fwdppConfig no longer attempts linkage to dependent libraries
* fixed errors in fwdpp/initms.hpp caught by clang++

## 0.4.0

* fwdpp/sugar/sampling.hpp added.  This streamlines taking samples from populations
* Unit test unit/sugar_sampling.cc added
* New test added to unit/crossoverTest.cc
* libsequence is no longer an installation dependency

## 0.3.9

This release fixes a minor bug for very high mutation rates.  Upgrading to this version is recommended!

* Fixed issues 26 and 27.
* New unit tests added to check sampling from populations

## 0.3.8

* Added some include directives whose omission resulted in failure to compile on clang.
* Doxygen input file now generated by ./configure.  No more need to manually update the version.

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
