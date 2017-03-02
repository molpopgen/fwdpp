## The fwdpp home page

This is the homepage for __fwdpp__, which is a C++ template library for implementing forward-time population-genetic simulations.

### What does fwdpp do?

__fwdpp__ provides the high-level abstractions for data types and algorithms that are needed to implement forward-time population-genetic simulations.  The abstractions are implemented using modern C++ template programming methods defined by the C++11 language standard.  

The goal of __fwdpp__ is to define fundamental _concepts_ in a generic way, such that different models may be implemented via minimal changes in code.  For example, data types have specific minimal API requirements, meaning that custom diploid, gamete, and/or mutation types may be easily implemented for your particular modeling needs.  Likewise, functions to generate new mutations, deterimine recombination breakpoints, etc., have simple minimum requirements.  In combination, custom data types + custom functions work together to generate specific models at compile time.  

The compile-time polymorphism described in the preceding paragraph can be further augmented using C++'s support for object-oriented programming.  By constructing custom base types for both types _and_ for functions (yes, C++ allows this!), one can also have run-time polymorphism, meaning that different models can be selected based on a command-line interface, for example.  

For more information, see the [overview](doc/md/overview) of library features.

#### fwdpy 

The combination of compile- and run-time polymorphism is used in our Python package, [fwdpy](http://molpopgen.github.io/fwdpy), which uses __fwdpp__ as its back-end.

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
