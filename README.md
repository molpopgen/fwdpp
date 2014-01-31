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

#Introduction

fwdpp is a C++ template library that abstracts the basic operations required to implement forward-time simulations of population- and quantitative-genetic models.  The library allows the simulation of single populations or metapopulations evolving under the standard evolutionary forces of drift, recombination, migration, and natural selection.  Arbitrary population size changes are also allowed. Different populations in a metapopulation may evolve under different fitness schemes.

The library uses advanced C++ techniques to allow arbitrary models to be implemented via the implementation of simple policies (see Documentation section below).  A programmer wishing to use the library will need a strong background in templates, function objects, and the Standard Template Library (STL).  Web resources for these topics vary too much in quality to recommend any particular one.  However, there are several classic books that are must-reads for C++ programmers (old school, I know):

1.  Scott Meyer's "trilogy" of "Effective C++", "More Effective C++", and "Effective STL".
2.  Nicolai Josuttis' "The C++ Standard Template Library"
3.  David Vandevoorde and and Nicolai Josuttis, "C++ Templates"

The first two are excellent books for people already familiar with C++ syntax but want to know more about effective software design using the language. Meyer's books are particularly good, espectially the first two.  The C++ Templates book is a bible of how to get the most out of templates.  It is a very advanced and detailed book, but I've found it helpful over the years.

The library user will also need some familiarity with the [boost](http://www.boost.org) libraries, especially "bind" and "function".  I refer the user to the boost website for the relevant documentation.

##Which C++?

fwdpp does not use any features from then newly-released C++11 standard.  The new standard extends/simplifies the language, and therefore I expect the current code base to be C++11-compliant. As compiler support for C++11 becomes more widespread, the library will likely start to use some of those features, which will drastically improve readability of some of the nastier bits of template wizardry.

##Citation

The manuscript describing fwdpp is currently on [arxiv](http://arxiv.org/abs/1401.3786).

#Dependencies

fwdpp depends upon the following libraries:

[boost](http://www.boost.org)<br>
[GSL](http://gnu.org/software/gsl)<br>
[zlib](http://zlib.net)<br>
[libsequence](http://github.com/molpopgen/libsequence)<br>

The first three are  available as pre-built packages on most Linux distributions.  The latter (libsequence) also depends on the first three, and must be built from source.

#Installation

./configure<br>
make<br>
make install<br>

##If dependent libraries are in non-stanard locations.

For example, if libsequence is in /opt:

CXXFLAGS=-I/opt/include LDFLAGS="$LDFLAGS -L/opt/lib" ./configure<br>
make<br>
make install

##Installing in a custom location

./configure --prefix=/path/2/where/you/want it

For example:

./configure --prefix=$HOME

#Documentation

The documentation is in the doc subdirectory that comes with the library.  There are two major pieces of documentation.  First is the detailed documentation of all library functions.  This is generated via [doxygen](http://www.doxygen.org), and the output is a folder called html.  To view the documentation, point a browser to html/index.html.  

The second piece of documentation is a tutorial on writing policies conforming to what fwdpp expects.  This document is doc/policies.tex and a pdf file of the documentation may be obtained by processing the file as follows:

cd doc<br>
pdflated policies<br>
pdflated policies<br>

One runs pdflatex twice to ensure that cross-references within the document are processed properly.