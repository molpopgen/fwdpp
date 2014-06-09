#Tests of the library performance

This folder contains various methods for testing that fwdpp gives correct output.

Many of these tests are meant to be run from within Rstudio using Rmarkdown (Rmd) files.  You will need R's knitr package installed for this to work.

Some of the longer tests will not use Rstudio, and realistically require a cluster to perform.

##Test 1: fixation probability

This test estimates the fixation probability of a mutation subject to genic selection.  It uses the example program pfix, which comes with the library.

To execute this test, process the file FixProb.Rmd from the folder fwdpp/test that comes with the library source code.  Be sure to compile the examples in fwdpp/examples first!

You can see the output that I get on my office machine by looking at the file [FixProbKRT.md](FixProbKRT.md).  (Note: if you run this example, you'll likely over-write the figure that my file is expecting to load!!)

This test compares an R implementation to a simple fwdpp program.  The best agreement to analytical predictions is for the smallest selection coefficients, as expected given the assumptions behind the analytical results.

##Test 2: compare neutral models to Hudson's ms.

This test is performed by executing NeutralModels.Rmd in fwdpp/test.

It assumes that you have [msstats](https://github.com/molpopgen/msstats) installed and in your user's path.  You also need Hudson's [ms](http://home.uchicago.edu:~/rhudson1) in your path.

The test plots the ECDF of four statistics calculated based on 1,000 replicates of ms and compares them to 100 replicates of the example program "diploid".  The parameters are theta = rho = 10.

The summary stats are:

1. S = number of segregating sites (Watterson 1975)
2. "Pi" or "sum of site heterogzygosity" = sum of 2pq over the S sites. (Nei, Tajima, others).
3. Minimum number of recombination events (Hudson and Kaplan)
4. H', a normalized version of Fay and Wu's H.

Two comparisons are made to diploid.  One is with N=100 and another with N = 1000.

This test merely eyeballs the agreement of the forward and backwards simulations.  Ideally, you'd do many more simulations and conduct statistical tests on the distributions (as was done in the fwdpp paper).  You can modify NSIMS and this test to do such analysis if you wish.

The output that I get (which also shows the seeds I used) is [here](NeutralModelsKRT.md).
