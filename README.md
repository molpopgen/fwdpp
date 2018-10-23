	fwdpp - A C++ template library for forward-time population genetic simulations



  Copyright (C) 2013 Kevin Thornton

  fwdpp is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Comments are welcome.

	- Kevin Thornton <krthornt@uci.edu>

# Preface

This README is the main page of the fwdpp documentation.  It may display some things oddly when viewed on GitHub and/or links to other documentation pages may be broken.  These issues are not bugs -- please see the Reference Manual at the [fwdpp home page](http://molpopgen.github.io/fwdpp/) for a fully-formatted version of this document.

# Using simulations in Python

If you are reading this, I'm guessing that you are more interested in _running_ forward simulations than you are in _developing_ them.  This package is about the latter.  I'd recommend that you check out [fwdpy11](http://molpopgen.github.io/fwdpy11), which is a Python package for running forward simulations and processing the output.  [fwdpy11](http://molpopgen.github.io/fwdpy11) is developed using fwdpp (this package).  [fwdpy11](http://molpopgen.github.io/fwdpy11) is obtainable from GitHub and will soon be avalable from [PyPi](http://pypi.python.org) or via the [bioconda](https://bioconda.github.io/) channel of [Anaconda](https://www.continuum.io/downloads).

fwdpy11 replaces [fwdpy](http://molpopgen.github.io/fwdpy), which was a Cython-based package for enabling fwdpp-based
simulations in Python.  We've end-of-life'd that project, but it is still on PyPi/Bioconda/GitHub for archival purposes,
as a couple of publications used it.  Trust us, though, [fwdpy11](http://molpopgen.github.io/fwdpy11) does everything
fwdpy did, but better.

# fwdpp-users mailing list

There is a [Google Group](https://groups.google.com/forum/#!forum/fwdpp-users) for fwdpp users.  It is called fwdpp-users.  I prefer that questions, etc., are raised there.  Questions on that list that are determined to be real issues will result in issue tickets on the GitHub repo for this site.  The group is also a good place for installation help, use questions, etc.

# Build status

* Status of master branch: [![Build Status](https://travis-ci.org/molpopgen/fwdpp.svg?branch=master)](https://travis-ci.org/molpopgen/fwdpp)
* Status of dev branch: [![Build Status](https://travis-ci.org/molpopgen/fwdpp.svg?branch=dev)](https://travis-ci.org/molpopgen/fwdpp)

# Introduction

fwdpp is a C++ template library that abstracts the basic operations required to implement forward-time simulations of population- and quantitative-genetic models.  The library allows the simulation of single populations or metapopulations evolving under the standard evolutionary forces of drift, recombination, migration, and natural selection.  Arbitrary population size changes are also allowed. Different populations in a metapopulation may evolve under different fitness schemes.

The library uses advanced C++ techniques to allow arbitrary models to be implemented via the implementation of simple policies (see Documentation section below).  A programmer wishing to use the library will need a strong background in templates, function objects, and the Standard Template Library (STL).  Web resources for these topics vary too much in quality to recommend any particular one.  However, there are several classic books that are must-reads for C++ programmers (old school, I know):

1.  Scott Meyer's "trilogy" of "Effective C++", "More Effective C++", and "Effective STL".
2.  Nicolai Josuttis' "The C++ Standard Template Library"
3.  David Vandevoorde and Nicolai Josuttis, "C++ Templates"

The first two are excellent books for people already familiar with C++ syntax but want to know more about effective software design using the language. Meyer's books are particularly good, espectially the first two.  The C++ Templates book is a bible of how to get the most out of templates.  It is a very advanced and detailed book, but I've found it helpful over the years.

## A note about which version to use

This code is distributed via my GitHub [account](http://www.github.com/molpopgen).  The "master" and "dev" branches should be viewed as experimental.  The [releases](https://github.com/molpopgen/fwdpp/releases), however, correspond to tested versions of the library fit for public consumption.  This means that, while the version number in the configure script on master/dev may match that of a recent release, _that does not mean that the features/stability/bugs present in master/dev are identical to those of the release._  If you want to use fwdpp for research, use the latest [release](https://github.com/molpopgen/fwdpp/releases).  If you want to play around with the latest and (occasionally not-so) greatest, look at the dev branch.  If you want to look at the latest I believe to be stable, look at master.  Also note that master may be ahead of dev, etc., depending on what I've committed from my development server to the repo stored at github.

### Revision history

[Release notes](@ref md_md_RELEASE_NOTES)

You should pay particular attention to any statements about backwards compatibility.  If you are updating from a previous version of __fwdpp__, you may wish to read this document regarding how to reproduce results based on previous versions: @ref md_md_preproc.

Specific version numbers ("tags" in git-ese, a.k.a. "releases") will occur when new feature are added to the library and/or bugs are fixed.  The details of what happens in each release can be found [here](@ref md_md_RELEASE_NOTES), beginning with release 0.2.4.


## Which C++?

As of version 0.2.5, fwdpp requires a compiler supporting the "C++11" version of the language.  Currently, fwdpp requires that your compiler support the flag -std=c++11 in order to use c++11 language features. Recent version of GCC  and clang both support this option, which covers most Linux and OS X users.

## Citation

The fwdpp manuscript has published in Genetics.  The accepted version of the manuscript is [here](http://www.genetics.org/content/early/2014/06/19/genetics.114.165019.abstract).  For LaTeX users:

~~~
@Article{,
  author = 	 {K. R. Thornton},
  title = 	 {A C++ Template Library for Efficient Forward-Time Population Genetic Simulation of Large Populations},
  journal = 	 {Genetics},
  year = 	 {2014},
  OPTkey = 	 {},
  volume = 	 {198},
  OPTnumber = 	 {},
  pages = 	 {157-166},
  OPTmonth = 	 {},
  OPTnote = 	 {},
  annote = 	 {doi:/10.1534/genetics.114.165019}
}
~~~

The version of fwdpp used in that publication is 0.2.4.

### Where fwdpp is now (compared to the publication).

The published version of fwdpp described a library with the following features:

* Objects (mutations and gametes) are stored only once in _doubly-linked lists_.
* Diploids are pairs of pointers to gametes (technically, pairs of iterators derived from the list of gametes).
* Gametes contain vectors of pointers (again, C++ iterators) to mutations.
* Each generation, extinct mutations and gametes are removed from the population.  This removal is a constant-time operation due to the use of linked lists.

This design had a certain elegance to it.  A pointer to a diploid gave you immediate access to the mutations by means of the pointer structure.  It was also compact in memory, because iterators are only 8 bytes (on a 64-bit system) while gametes and mutations are 64 and _at least_ 48 bytes, respectively.

The fwdpp publication showed that the library performs well in terms of speed compared to other tools out there.  However, the original design had the following problems:

* Several lookup operations were implemented using expensive linear-time searches.
* The linked lists plus the constant insertion and deletion of new and extinct objects, respectively, resulted in memory fragmentation or poor "cache locality".

The 0.3.x releases of fwdpp solved most of the first problem, and sped the library up by over an order of magnitude.  To address the second problem (poor cache locality), I recommended that users link programs based on fwdpp to an external library replacing the built-in malloc (the main memory allocation function for the C family of languages, and C++'s "new" is a wrapper around malloc).  I specifically recommended using Google's tcmalloc.  It turned out that simply using an industrial-strength memory allocator went a long way towards addressing the performance hit due to memory fragmentation.

(These releases also introduced various sub-libraries aimed at making fwdpp easier to use.)

At this point, intuition (backed up by extensive code profiling using tools like Google's profiling library and valgrind's cachegrind) told me that the main hurdle to improving performance was to address memory usage.

Release 0.4.4 introduced a fundamental change in the library design:

* Extinct objects are no longer deleted.  Instead, their locations are recorded in a FIFO queue.  These queues are used to "recycle" already-allocated memory.
* (Similarly, fixed variants can be tagged for recycling under many modeling situations.)
* Once extinct objects are no longer removed, then there are no longer insertions/removals to/from the middle of containers.  Rather, we either recycle an object, or add a newly-allocated object to the end.
* Thus, we can replace doubly-linked lists with vectors.
* Further, we can replace iterators with integers.

All of these changes were introduced in one fell swoop in 0.4.4, along with a set of other API changes that I'd wanted to make for a while.  The current design is conceptually the same as the published version, with 8 byte keys representing where mutations and gametes are, thus ensuring that an object is only represented once in memory.

However, the new design is also less elegant.  Now, the vectors of gametes and mutations have to be passed along with the diploids.

So, why do this? __It is a lot faster!__  Simulations of large genomic regions in large populations can be up to 80% faster!  In fact, tcmalloc isn't necessary to get really good performance any more.  Using it still improves run-times by about 10% (on Intel systems at least...), but that isn't a lot compared to the 50% improvement that it gave to previous versions of the library.


# Documentation

## Online

A tutorial on policies and the library's reference manual can be found at [molpopgen.github.io/fwdpp](http://molpopgen.github.io/fwdpp) or [here](@ref md_md_policies).

__Note:__ the links above may be out of date, as the online documentation are not regenerated automatically.  If you want the latest, builds the docs from source.

## Tutorials

__The first few tutorials are out of date due to API changes in 0.4.4.  The basic ideas are the same, but you'll have to look at the example programs and unit tests for the best HOW-TOs until I fix these documents.__

The [fwdpp](http://molpopgen.github.io/fwdpp) main page contains several tutorials:

* @ref md_md_datatypes
* @ref md_md_policies
* @ref md_md_multiloc
* @ref md_md_serialization
* @ref md_md_devtools
* @ref md_md_sugar


## Built from source
The source code documentation is in the doc subdirectory that comes with the library.  There are two major pieces of documentation.  First is the detailed documentation of all library functions.  This is generated via [doxygen](http://www.doxygen.org), and the output is a folder called html.  To view the documentation, point a browser to html/index.html.

## Example documentation
The examples can be read in html form via the online reference manual linked to above.  You can find the two simplest examples online at the fwdpp [wiki](https://github.com/molpopgen/fwdpp/wiki) on github.

# Projects using fwdpp:

* [fwdpp_perf](https://github.com/molpopgen/fwdpp_perf) is a collection of programs showing how to run independent simulations in a multi-core/many-core environment.  Example programs use either C++11 threads or MPI to run simulations.  I also use this package for performance testing/benchmarking/profiling/etc.
* [fwdpy11](https://github.com/molpopgen/fwdpy11) brings fwdpp-powered simulations to the Python programming language.


# Dependencies

## System requirements

You must have the following on your system:

1. A C++ compiler that supports the C++11 language standard.
2. Ideally, one should have the [git](http://git-scm.com/book/en/Getting-Started-Installing-Git) command line tools installed.  These are likely already installed on many systems.

We routinely test fwdpp on several different systems with different compilers.  The Travis CI build matrix for fwdpp includes:

* A miniconda3 environment using gcc 4.8.5 with [libsequence](http://github.com/molpopgen/libsequence) installed via [bioconda](https://bioconda.github.io/).
* GCC 5 or 6 on Ubuntu, with either C++11 or C++14.  [libsequence](http://github.com/molpopgen/libsequence) is not used for these builds.

## Library dependencies

The minimal dependencies required to use the library to develop simulations are:

1.  [GSL](http://gnu.org/software/gsl)
2.  [zlib](http://zlib.net)

In order to compile the example programs, the following dependency is optional as of fwdpp 0.5.7: 

1.  [libsequence](http://github.com/molpopgen/libsequence).

As of fwdpp 0.7.0, some example programs require the `program_options` library from [boost](http://www.boost.org).  This
is one of boost's runtime libraries, meaning that you need a "full" boost installation. Conda provides such a thing, and
Linux variants like Ubuntu provide the different runtime libraries as separate packages.

The configure script checks for these optional dependencies and will skip compilation of various example programs if
they are missing.

In order to compile the unit tests, you also need:

1.  [boost](http://www.boost.org).

For OS X users, all of the above dependencies are available via [homebrew](http://brew.sh) or [conda](http://conda.pydata.org/docs/)/[bioconda](https://bioconda.github.io).

## Performance

Performance testing has been moved to the [fwdpp_perf](http://github.com/molpopgen/fwdpp_perf) project.

## Obtaining the source code

### Obtaining the master branch
You have a few options:

1. Clone the repo (best option): git clone https://github.com/molpopgen/fwdpp.git
2. Click on "Download Zip" at https://github.com/molpopgen/fwdpp


### Obtaining a specific release
Again, a few options:

1. Click on "Releases" at https://github.com/molpopgen/fwdpp, then download the one you want
2.  Clone the repo (see previous section)

Or:

1. Get a list of releases by saying "git tag -l"
2. Checkout the release you want.  For example "git checkout 0.2.6"

Another option is to click on [releases](https://github.com/molpopgen/fwdpp/releases) from the main [project page](https://github.com/molpopgen/fwdpp) at github.

On a decent browser, when you click on a release, it should be called fwdpp-version.tar.gz.  Sometimes, though, you may get version.tar.gz.  This is a browser-by-github interaction problem.  On my systems, I get the correct result.

# Installation

## What does fwdpp install?

Two things:

* The library itself, which is simply a set of C++ header files
* A single binary called fwdppConfig, which you may use to check what version you have installed on your system.

~~~{sh}
fwdppConfig --version
~~~

fwdppConfig was introduced in version 0.3.3.  Its main raison d'etre is in helping other configure scripts test for what version is on their system.  For example, which may go into a configure.ac file for your project:

~~~{sh}
dnl Check that fwdpp version is sufficient
if test "$FWDPPVERSION" \> "0.3.2"
then
	echo "fwdpp version $FWDPPVERSION detected"
else
	AC_MSG_ERROR([fwdpp >= 0.3.3 required, please install from http://github.com/molpopgen/fwdpp])
fi
~~~

## The case of a standard system with all dependencies installed in standard locations

If you cloned the git repo:
~~~
cd fwdpp
~~~
If you downloaded a release:

~~~
tar xzf fwdpp-version.tar.gz
cd fwdpp-version
~~~

Then:

~~~
./configure
make
make install
~~~

To compile examples and unit tests, and to execute unit tests:

~~~
make check
~~~

Currently, the example programs will not get installed via "make install".   If you want them to be installed system-wide, copy the binaries manually to where you need them.

## To compile examples and unit tests

__Note:__ if you only wish to compile the example programgs, issue the 'make check' command from the example subdirectory.  This will allow users without boost on their system to compile the examples but not attempt to compile the unit tests (which will fail to compile on systems without boost).

You will need [libsequence](http://github.com/molpopgen/libsequence) installed on your system in order to compile the example programs.

## If dependent libraries are in non-stanard locations.

For example, if libsequence is in /opt:

~~~{.sh}
#Note, you need to add in the desired optimization (-OXX) level:
CPPFLAGS="-I/opt/include" LDFLAGS="$LDFLAGS -L/opt/lib -Wl,-rpath/opt/lib" ./configure 
make check
make install
~~~

The "-Wl,-rpath" business adds the library path locations into the compiled binaries.  The flags work on GCC and clang
compilers, and are often handy when using systems like Anaconda or Homebrew to manage dependencies.

## Installing in a custom location

~~~
./configure --prefix=/path/2/where/you/want it
~~~
For example:

~~~
./configure --prefix=$HOME
~~~

## Examples

The examples are documented [here](@ref md_md_examples)
