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