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

# Introduction

fwdpp is a C++ template library that abstracts the basic operations required to implement forward-time simulations of population- and quantitative-genetic models.  The library allows the simulation of single populations or metapopulations evolving under the standard evolutionary forces of drift, recombination, migration, and natural selection.  Arbitrary population size changes are also allowed. Different populations in a metapopulation may evolve under different fitness schemes.

The library uses advanced C++ techniques to allow arbitrary models to be implemented via the implementation of simple policies (see Documentation section below).  A programmer wishing to use the library will need a strong background in templates, function objects, and the Standard Template Library (STL).  Web resources for these topics vary too much in quality to recommend any particular one.  However, there are several classic books that are must-reads for C++ programmers (old school, I know):

1.  Scott Meyer's "trilogy" of "Effective C++", "More Effective C++", and "Effective STL".
2.  Nicolai Josuttis' "The C++ Standard Template Library"
3.  David Vandevoorde and and Nicolai Josuttis, "C++ Templates"

The first two are excellent books for people already familiar with C++ syntax but want to know more about effective software design using the language. Meyer's books are particularly good, espectially the first two.  The C++ Templates book is a bible of how to get the most out of templates.  It is a very advanced and detailed book, but I've found it helpful over the years.

##A note about which version to use

This code is distributed via my gitub [account](http://www.github.com/molpopgen).  The "master" and "dev" branches should be viewed as experimental.  The [releases](https://github.com/molpopgen/fwdpp/releases), however, correspond to tested versions of the library fit for public consumption.  This means that, while the version number in the configure script on master/dev may match that of a recent release, _that does not mean that the features/stability/bugs present in master/dev are identical to those of the release._  If you want to use fwdpp for research, use the latest [release](https://github.com/molpopgen/fwdpp/releases).  If you want to play around with the latest and (occasionally not-so) greatest, look at the dev branch.  If you want to look at the latest I believe to be stable, look at master.  Also note that master may be ahead of dev, etc., depending on what I've committed from my development server to the repo stored at github.

###Revision history

[Release notes](@ref md_md_RELEASE_NOTES)

You should pay particular attention to any statements about backwards compatibility.  If you are updating from a previous version of __fwdpp__, you may wish to read this document regarding how to reproduce results based on previous versions: @ref md_md_preproc.

Specific version numbers ("tags" in git-ese, a.k.a. "releases") will occur when new feature are added to the library and/or bugs are fixed.  The details of what happens in each release can be found [here](@ref md_md_RELEASE_NOTES), beginning with release 0.2.4.

##Which C++?

As of version 0.2.5, fwdpp requires a compiler supporting the "C++11" version of the language.  Currently, fwdpp requires that your compiler support the flag -std=c++11 in order to use c++11 language features. Recent version of GCC (4.7 or greater) and clang (3.4 or greater, but I've not checked earlier versions) both support this option, which covers most Linux and OS X users.

##Citation

The fwdpp manuscript has been accepted for publication in Genetics.  The accepted version of the manuscript is [here](http://www.genetics.org/content/early/2014/06/19/genetics.114.165019.abstract).  For LaTeX users:

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


# Documentation

##Online

A tutorial on policies and the library's reference manual can be found at [molpopgen.github.io/fwdpp](http://molpopgen.github.io/fwdpp) or [here](@ref md_md_policies).

__Note:__ the links above may be out of date, as the online documentation are not regenerated automatically.  If you want the latest, builds the docs from source.

## Tutorials

The [fwdpp](http://molpopgen.github.io/fwdpp) main page contains several tutorials:

* @ref md_md_datatypes
* @ref md_md_policies
* @ref md_md_multiloc
* @ref md_md_serialization
* @ref md_md_devtools
* @ref md_md_sugar


##Built from source
The source code documentation is in the doc subdirectory that comes with the library.  There are two major pieces of documentation.  First is the detailed documentation of all library functions.  This is generated via [doxygen](http://www.doxygen.org), and the output is a folder called html.  To view the documentation, point a browser to html/index.html. 

##Example documentation
The examples can be read in html form via the online reference manual linked to above.  You can find the two simplest examples online at the fwdpp [wiki](https://github.com/molpopgen/fwdpp/wiki) on github.

# Dependencies

##System requirements

You must have the following on your system:

1. A C++ compiler equivalent to gcc 4.6 or greater.  On OS X, this likely means that you should be running OS X Mavericks with Xcode installed, and have then installed the command line tools (which is done from within Xcode.)
2. Ideally, one should have the [git](http://git-scm.com/book/en/Getting-Started-Installing-Git) command line tools installed.  These are likely already installed on many systems.

##Library dependencies

The minimal dependencies required to use the library to develop simulations are:

1.  [GSL](http://gnu.org/software/gsl) 
2.  [zlib](http://zlib.net) 
3.  [libsequence](http://github.com/molpopgen/libsequence).

In order to compile the unit tests, you also need:

1.  [boost](http://www.boost.org).

In order to get good performance out of simulations you also need one or both of the following:

1.  [boost](http://www.boost.org).
2.  [Google's perftools](https://code.google.com/p/gperftools/)

For OS X users, all of the above dependencies are available via [homebrew](http://brew.sh).

### Why use boost and/or Google's perftools?

Forward simulations are constantly allocating relativel small objects whose lifetimes are random variables.  For example, mutations enter populations every generation and are typically rapidly lost.   We can improve the performance of a simulation using various techniques, including replacing the C++ standard library's allocator with a custom allocator and/or replacing the C library's malloc function with a faster replacement.

The boost library provides fast_pool_allocator, which is a good replacement for the standard allocators.  Google's perftools provide the tcmalloc library, which provides a faster malloc.

Empirically, on my linux systems, I have found that the fastest combination appears to be the standard C++ allocator combined with linking to tcmalloc.  See below for how to compile the examples in various ways, allowing you to benchmark performance on your own system.

[See here](performance/tcmalloc.pdf) for an example of how tcmalloc can improve performance.  Results are based on running the example program diploid_ind with N=1000, theta=rho=5000 for 10N generations.

##Obtaining the source code

###Obtaining the master branch
You have a few options:

1. Clone the repo (best option): git clone https://github.com/molpopgen/fwdpp.git
2.  Click on "Download Zip" at https://github.com/molpopgen/fwdpp 


###Obtaining a specific release
Again, a few options:

1. Click on "Releases" at https://github.com/molpopgen/fwdpp, then download the one you want
2.  Clone the repo (see previous section)

Or:

1. Get a list of releases by saying "git tag -l"
2. Checkout the release you want.  For example "git checkout 0.2.6"

Another option is to click on [releases](https://github.com/molpopgen/fwdpp/releases) from the main [project page](https://github.com/molpopgen/fwdpp) at github.

On a decent browser, when you click on a release, it should be called fwdpp-version.tar.gz.  Sometimes, though, you may get version.tar.gz.  This is a browser-by-github interaction problem.  On my systems, I get the correct result.

# Installation

##The case of a standard system with all dependencies installed in standard locations

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

To compile examples and unit tests:

~~~
make check
~~~

Currently, the example programs will not get installed via "make install".   If you want them to be installed system-wide, copy the binaries manually to where you need them.

##To compile examples and unit tests

__Note:__ if you only wish to compile the example, issue the 'make check' command from the example subdirectory.  This will allow users without boost on their system to compile the examples but not attempt to compile the unit tests (which will fail to compile on systems without boost).

### Using boost containers and the standard C-library malloc

This is the default:

~~~{.sh}
./configure
make check
~~~

### Using the standard C++ containers and the standard C library malloc

~~~{.sh}
./configure --enable-standard=yes
make check
~~~

### Using boost containers and Google's tcmalloc replacement for malloc

~~~{.sh}
./configure --enable-tcmalloc=yes
make check
~~~


### Using the standard C++ library containers and Google's tcmalloc replacement for malloc

~~~{.sh}
./configure --enable-standard=yes --enable-tcmalloc=yes
make check
~~~~

Note:  none of these options affect any other programs that you write using fwdpp!  You'll have to manage that on your own for your projects' builds.

##If dependent libraries are in non-stanard locations.

For example, if libsequence is in /opt:

~~~{.sh}
#Note, you need to add in the desired optimization (-OXX) level:
./configure CXXFLAGS=-"-O2 -I/opt/include" LDFLAGS="$LDFLAGS -L/opt/lib"
make check
make install
~~~

##Installing in a custom location

~~~
./configure --prefix=/path/2/where/you/want it
~~~
For example:

~~~
./configure --prefix=$HOME
~~~

##Examples

The examples are documented [here](@ref md_md_examples)
