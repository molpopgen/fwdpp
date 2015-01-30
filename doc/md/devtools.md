# Tutorial 4: Starting projects based on fwdpp

## Intro

This document covers how to use GNU autotools to manage a project based on [fwdpp](http://molpopgen.github.io/fwdpp/).

## What you need

You'll need the following installed on your system

* [autoconf](https://www.gnu.org/software/automake/)
* [automake](https://www.gnu.org/software/automake/)

For Linux users, these are probably already installed on your system.

## Getting the right tools on OS X

OS X users are at a slight disadvantage here, as the Xcode tool suite no longer includes the GNU tools setup.  However, you can use [homebrew](http://brew.sh) to get these tools.  Please see that website for complete instructions, but at the time of this writing, the following will work:

1. Make sure that you have Xcode installed.  You cannot in stall brew without it.
2. Install brew:

~~~{.sh}
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
~~~

3. Install several packages that you'll need:

~~~{.sh}
brew install autoconf automake libtool wget
~~~

4. Optionally, install some other handy tools:

~~~{.sh}
brew install doxygen
~~~

## Generating a skeleton project

Move into the devtools directory of the fwdpp source code repository:

~~~{.sh}
cd devtools
~~~

There, you'll find an executable bash-language script called setup.sh.  It's usage is:

~~~{.sh}
Usage: ./setup.sh options
Mandatory options are:
 -d/--dir root path to place the project (DIR)
 -p/--project Name for project directory (PROJECT)
Together, the -p and -d options will create a skeleton package in DIR/PROJECT.
For example: ./setup.sh -p ~/src -d my_fwdpp_project
Optional options are:
 -u/--url project url.  This must be single-quoted with special characters escaped, e.g. -u 'https:\/\/github.com\/molpopgen\/fwdpp'
~~~

Here is an example of using this script to create a project called fwdpp_project in my user's src subdirectory.  I intend for the project to be hosted on [GitHub](http://github.com), so I'll pass it an appropriate URL, too:

~~~{.sh}
./setup.sh -d ~/src -p fwdpp_project -u 'https:\/\/github.com\/molpopgen\/fwdpp_project'
~~~

The following got printed to the screen:

~~~
Base path set to '/Users/krthornt/src'
Project name set to 'fwdpp_project'
~~~

We're almost done, but our build system isn't 100 percent ready for our system.  To finish the job, we need to say:

~~~{.sh}
cd ~/src/fwdpp_project
autoreconf -fi
autoheader
automake --add-missing -c
~~~

Hint: at this point, it may be prudent to make an initial commit to your version control system!

## Some details

On all systems that I work on, the build system now works after the commands in the previous section are completed.

Our configure.ac contains the correct project name, URL, etc:

~~~
head -n 10 configure.ac 
AC_PREREQ(2.59)

AC_INIT([fwdpp_project],[0.1.0],[https://github.com/molpopgen/fwdpp_project])
AC_CONFIG_SRCDIR([src/fwdpp_project.cc])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
~~~

Further, src/Makefile.am looks good:

~~~
bin_PROGRAMS=fwdpp_project

fwdpp_project_SOURCES=fwdpp_project.cc

AM_CPPFLAGS=-W -Wall -I. -I..

if DEBUG
else !DEBUG
AM_CPPFLAGS+=-DNDEBUG
endif
~~~

The configure script checks for the presence of several [boost](http://www.boost.org) libraries, [fwdpp](http://molpopgen.github.io/fwdpp/), [libsequence](http://molpopgen.github.io/libsequence), the [GSL](http://www.gnu.org/software/gsl/), and [zlib](http://zlib.net).  It also checks that your compiler is C++11-ready.

## What next?

Now, you are basically on your own, and you need to implement your simulation, etc.  See the various other pieces of documentation:

* @ref md_md_policies
* @ref md_md_multiloc
* @ref md_md_serialization

Also, take a look at the source code for the example programs that come with the library
