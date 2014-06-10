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

The test plots the ECDF of four statistics calculated based on 1,000 replicates of ms and compares them to 100 replicates of the example program "diploid_ind".  The parameters are theta = rho = 10.

The summary stats are:

1. S = number of segregating sites ([Watterson, 1975](http://www.ncbi.nlm.nih.gov/pubmed/1145509))
2. "Pi" or "sum of site heterogzygosity" = sum of 2pq over the S sites. (Nei, [Tajima](http://www.genetics.org/content/105/2/437.abstract), others).
3. Minimum number of recombination events ([Hudson and Kaplan](http://www.genetics.org/content/111/1/147.abstract))
4. [Zeng et al.'s H'](http://www.genetics.org/content/174/3/1431.abstract), a normalized version of [Fay and Wu's H](http://www.genetics.org/content/155/3/1405.abstract).

Two comparisons are made to diploid_ind.  One is with N=100 and another with N = 1000.

The output that I get (which also shows the seeds I used) is [here](NeutralModelsKRT.md).

##Test 3: a more rigorous comparison of neutral models to Hudson's ms.

This test is executed via NeutralModelsBetter.Rmd.  It is the same idea as the preceeding, except that 1,000 forward simulations are run and the p-values of Kolmogorov-Smirnov tests comparing the distributions of summary statistics between fwdpp and ms are printed.  Additionally, N=10,000 is also simulated.

__Warning:__ this test will take hours-to-days to run.  The N=10,000 case really should be done on a cluster.

Typically, you will see more differences between fwdpp and ms for smaller N.  This is expected.  Recall that ms simulates a sample of size n taken from an infinitely-large Wright-Fisher population.  In contrast, the forward simulation is taking a sample of the same size n from a much smaller population.  However, as the N in the forward simulation increases, one should see convergence in outcomes to the predictions of the infinite-N models.

It is very important to understand that the 21 p-values reported for each value of N are __not__ the outcomes of indepdent statistical tests.   For example, statisics 1,3, and 4 in the list in the preceeding section are highly correlated with one another, and correlations exists between any pair of statistics summarizing variation data.  Thus, the rate at which the K-S test will reject the null model cannot be assumed to be the usual alpha.

###Running this test on a cluster

If you have a Gride Engine (GE) queuing system, you can run the fwdpp part of this test using a script like the following.  (Note: this script is bad practice, as it writes each replicate to a separate file.  It would be much smarter to pipe the output through a [file locking tool](https://github.com/molpopgen/atomic_locker), but I skip that here in the interest of simplicity.)

These scripts assume that all binaries are available in your user's path, and that your cluster is setup to automatically cd to SGE_O_WORKDIR.

You will need to replace "queuename" with the appropriate queue names for your cluster, add lines to add relevant modules, etc.

Job 1 runs the simulations:

```{sh}
#!/bin/sh

#$ -q queuename
# Set up array job of 1,000 job
#$ -t 1-1000
#Name the job
#$ -N job1
#Define model params
N=1000
G=`echo "10*$N"|bc`
#give each task unique seed.  This method controls for 2 jobs starting at same time such that RANDOM returns same seed
SEED=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`

diploid_ind $N 10 10 $G 10 1 $SEED | gzip > out.$N.$SGE_TASK_ID.gz
```

Job 2 cleans up after job 1:

```{sh}
#$/bin/sh

#$ -N cleanup
#Do not execute in queue until job 1 is done
#$ -hold_jid job1

N=1000
zcat out.$N.*.gz | msstats | gzip > msstats.$N.gz
```

You can them process the output in R as done in the Rmd file for this test.

Users of torque/pbs will need to figure out how to do the above on their own.
