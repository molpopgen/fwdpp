#Tests of the library performance

This folder contains various methods for testing that fwdpp gives correct output.

Many of these tests are meant to be run from within Rstudio using Rmarkdown (Rmd) files.  You will need R's knitr package installed for this to work.

Some of the longer tests will not use Rstudio, and realistically require a cluster to perform.

##Test 1: fixation probability

This test estimates the fixation probability of a mutation subject to genic selection.  It uses the example program pfix, which comes with the library.

To execute this test, process the file FixProb.Rmd from the folder fwdpp/test that comes with the library source code.  Be sure to compile the examples in fwdpp/examples first!

You can see the output that I get on my office machine by looking at the file [FixProbKRT.md](FixProbKRT.md).  (Note: if you run this example, you'll likely over-write the figure that my file is expecting to load!!)
