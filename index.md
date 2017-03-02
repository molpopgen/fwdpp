## The fwdpp home page

This is the homepage for __fwdpp__, which is a C++ template library for implementing forward-time population-genetic simulations.

### What does fwdpp do?

__fwdpp__ provides the high-level abstractions for data types and algorithms that are needed to implement forward-time population-genetic simulations.  The abstractions are implemented using modern C++ template programming methods defined by the C++11 language standard.  

The goal of __fwdpp__ is to define fundamental _concepts_ in a generic way, such that different models may be implemented via minimal changes in code.  For example, data types have specific minimal API requirements, meaning that custom diploid, gamete, and/or mutation types may be easily implemented for your particular modeling needs.  Likewise, functions to generate new mutations, deterimine recombination breakpoints, etc., have simple minimum requirements.  In combination, custom data types + custom functions work together to generate specific models at compile time.  

The compile-time polymorphism described in the preceding paragraph can be further augmented using C++'s support for object-oriented programming.  By constructing custom base types for both types _and_ for functions (yes, C++ allows this!), one can also have run-time polymorphism, meaning that different models can be selected based on a command-line interface, for example.  

For more information, see the [overview](doc/md/overview) of library features.

#### fwdpy 

The combination of compile- and run-time polymorphism is used in our Python package, [fwdpy](http://molpopgen.github.io/fwdpy), which uses __fwdpp__ as its back-end.

### Using fwdpp to develop software.

You can use __fwdpp__ to develop several different types of applications.  The simplest type would be a command-line application.  The examples in the GitHub repo for this project are all command-line applications.  However, complex command-line programs quickly get unwieldy to implement and also to use.  It is useful to be able to run simulations in an interpreted environment such as Python or R.  Using Python, there are several possible options:

* [boost.python](http://www.boost.org/doc/libs/1_63_0/libs/python/doc/html/index.html).  This is one of the oldest means of exposing C++ types to Python.  The python_examples directory of the __fwdpp__ GitHub repo contains two examples.
* [pybind11](http://pybind11.readthedocs.io/) is a relatively new project attempting to provide a modern (C++11) library to get C++ and Python to work together.  It is a very nice project, and essentially replaces boost.python. 
* [Cython](http://www.cython.org) is a static compiler for Python that allows C++ types to be exposed.  This is what we use to develop [fwdpy](http://molpopgen.github.io/fwdpy), which exposes a lot of __fwdpp__ to Cython.  Those exposed types can be re-used in other projects, or [fwdpy](http://molpopgen.github.io/fwdpy) itself can be used as a base for your project.

The main difference between boost.python/pybind11 and Cython is that the former provide a "C++-first" perspective on writing Python extension modules.  That is, you can write your Python module in idiomatic C++.  However, that does not provide the most "Pythonic" experience for the folks using your package.  Cython, on the other hand, is Python-first, which causes some frustrations when wrapping C++ code, but ultimately leads to a more extensible package and a Pythonic experience.

It is also possible to use __fwdpp__ from [R](http://r-project.org) via the excellent [Rcpp](http://www.rcpp.org) package.  In fact, [fwdpy](http://molpopgen.github.io/fwdpy) started out as an R package, but we migrated to Cython/Python because Cython's grammar makes it easier for users to write their own custom extensions.



### Citation

If you use __fwdpp__ for yor research, either to develop simulations or you use the example programs, please cite the following manuscript:

* Thornton, K. R. (2014) A C++ template library for efficient forward-time population genetic simulation of large populations.  Genetics 98:157-166  PMID: 24950894, [Manuscript](http://www.genetics.org/content/198/1/157.abstract), [Software](https://github.com/molpopgen/fwdpp)

Here it is in Bibtex format:

~~~
@ARTICLE{Thornton2014-hx,
  title       = "A C++ template library for efficient forward-time population
                 genetic simulation of large populations",
  author      = "Thornton, Kevin R",
  affiliation = "Department of Ecology and Evolutionary Biology, University of
                 California, Irvine, California 92697",
  abstract    = "fwdpp is a C++ library of routines intended to facilitate the
                 development of forward-time simulations under arbitrary
                 mutation and fitness models. The library design provides a
                 combination of speed, low memory overhead, and modeling
                 flexibility not currently available from other forward
                 simulation tools. The library is particularly useful when the
                 simulation of large populations is required, as programs
                 implemented using the library are much more efficient than
                 other available forward simulation programs.",
  journal     = "Genetics",
  volume      =  198,
  number      =  1,
  pages       = "157--166",
  month       =  sep,
  year        =  2014,
  keywords    = "population genetics; quantitative genetics; simulation",
  language    = "en"
}
~~~

### Applications of fwdpp 

The library has been used for the following software projects:

* [fwdpy](https://github.com/molpopgen/fwdpy) uses __fwdpp__ to provide an environment for foward simulation in Python.  This is a very large-scale project, and is where a lot of our future forward-simulation-related tools will end up.

### Publications (that we are aware of...) using fwdpp

* Sanjak, Jaleal S., Anthony D. Long, and Kevin R. Thornton. 2017. “A Model of Compound Heterozygous, Loss-of-Function Alleles Is Broadly Consistent with Observations from Complex-Disease GWAS Datasets.” PLoS Genetics 13 (1): e1006573. [Paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006573)  This paper also uses [fwdpy](http://github.com/molpopgen/fwdpy) for some of the supplementary material.
* Beissinger, Timothy M., Li Wang, Kate Crosby, Arun Durvasula, Matthew B. Hufford, and Jeffrey Ross-Ibarra. 2016. “Recent Demography Drives Changes in Linked Selection across the Maize Genome.” Nature Plants 2 (June): 16084. [Paper](http://www.nature.com/articles/nplants201684?WT.feed_name=subjects_next-generation-sequencing)

### Early applications of fwdpp

Early version of the code that eventually became __fwdpp__ were used in the following publications:

* [Thornton et al.](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003258) on modeling "rare alleles of large effect" and genome-wide association studies
* [Baldwin-Brown et al.](http://mbe.oxfordjournals.org/content/31/4/1040.full) on modeling "evolve and resequence" experiments"
