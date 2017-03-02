# Example programs

The best documentation of how to use the library are the example simulations in the examples/ subdirectory of the source code repo.

###Running the examples:

This following list of parameters and their definitions is in common to all of the example programs:

* N = the number of diploid individuals to simulate
* theta = 4Nu, where u is the mutation rate per gamete per generation. For example, theta = 100 means that on average, 2Nu = 50 new mutations are entering the population each generation.
* rho = 4Nr, where r is the recombination rate per diploid per generation
* g = the number of generations to simulate. Often, this should be >= 8N at a minimum
* n = the sample size to draw at the end of the simulation. To match the typical modeling assumptions of population genetics, you should have n << 2N.
* nreps = the number of replicates to simulate
* seed = a random number seed. I will use $RANDOM as a seed, referring to the bash shell method to return a random integer.

Note: familiarity with Hudson's "[ms](http://home.uchicago.edu/~rhudson1)" program is helpful for some of what comes below.

####diploid_ind (diploid_ind.cc)

The first program, diploid, simulates a Wright-Fisher population with mutation, recombination, and drift. To run it:

~~~
./diploid_ind N theta rho g n nreps seed | gzip > outfile.gz
~~~

The data in outfile.gz will be in the same format as Dick Hudson's "ms" coalescent simulator.

Example:

~~~
./diploid_ind 10000 10 10 100000 50 1 $RANDOM | gzip > test_diploid.gz
~~~


####diploid\_binaryIO_ind (diploid_binaryIO_ind.cc)

The next program is called diploid\_binaryIO\_ind. This program is identical to diploid, except that it only simulates one replicate at a time, creates two output files, and outputs the entire population rather than a sample of size n << 2N. The first file is an index file, containing an integer representing the replicate number of the output and the position in the haplotypes file where this record begins. The second file is the haplotypes file, which contains the entire population in binary format. 

Usage:

~~~
./diploid_binaryIO N theta rho g replicate\_number index\_filename haplotype\_filename seed.
~~~

Example for Open Grid Engine compute clusters:

~~~{.sh}
#!sh
#$ -t 1-100
#$ -N DIPBIN

seed=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`
~~~

\#note: the below assumes that the binary is in the users's $PATH and that the GE system knows how to link it to the relevant dynamic libraries

~~~
diploid_binaryIO_ind 10000 10 10 100000 $SGE_TASK_ID indexfile hapfile $seed
~~~

The above script, when submitted to a Grid Engine queue, will result in 100 populations of size N=10,000 being written to hapfile. Further, “indexfile” will contain the ID number and position of each file. Records are not over-written because the program uses POSIX file locking to ensure that only 1 process at a time can do the writing. This is a complex program, as it mixes C++ objects with output streams such that they can be written to C-style file descriptors, which is required in order to use file locking (which is a C feature with no C++ analog). However, the advantage is that you write all data to one large file, avoiding the plague of lots of small files that can bring distributed file systems to their knees.


####diploid\_fixed\_sh\_ind (diploid_fixed_sh\_ind.cc) and diploid\_fixed\_sh\_ind_lambda (diploid_fixed_sh_ind_lambda.cc)

This program is similar to diploid_ind, but adds an additional mutation rate (theta\_selected = 4Nu\_s, where u\_s is the mutation rate per gamete per generation to mutations with selection coefficient s) to mutations with selection coefficient s and dominance h. Fitness across sites is multiplicative. The output is in "ms" format--one block for neutral mutations followed by one block for selected mutations.

Usage:

~~~
./dipoid_fixed_sh_ind N theta theta_selected rho s h g n nreps seed | gzip > outfile.gz
~~~

For this program, s can be positive or negative, as can h.

####diploid\_twopop\_mig\_ind (diploid_twopop_mig_ind.cc)

This program simulates an ancestral population for g generations, at which point a daughter population of size N is “budded” off from the ancestral population. Evolution continues for g2 more generations, with symmetric migration at rate M = 4Nm, where m is the migration rate per diploid per generation. The output is in "ms" format, with n haplotypes per sample just like how ms outputs data for multi-population models.

Usage:

~~~
./diploid_twopop_mig_ind N theta rho g g2 M n nreps seed | gzip > outfile.gz
~~~

Note: the demographic model here implemented may be viewed as biologically bizarre, as it mimics the default behavior of “ms” for population split models. Let’s compare a specific example vs. the equivalent ms command line:

~~~
./diploid_twopop_mig 10000 50 50 100000 1000 1 50 1 $RANDOM | gzip > outfile.gz
~~~

and

~~~
ms 100 1 -t 50 -r 50 1000 -I 2 2 1 -ej 0.025 2 1 -em 0.025 1 2 0.
~~~

Why may this be considered odd? In ms, when two populations are merged, the rate of coalescence is unaffected by default (this behavior is documented, and it is up to the user to adjust population sizes when population merge and split in ms). That means when the two populations, each of size N merge, the merged (ancestral) population is still of size N. diploid\_twopop\_mig is doing the same thing forwards in time: an ancestral population of size N magically changes into two populations of size N. 

####migsel\_ind (migsel_ind.cc)

Simulates 2 equal-sized populations of size N (N remains constant over time) diploids with selection at strength s and dominance h. Migration occurs between the two populations.

_Note that the two populations never share a common ancestor in this simulation!  In other words, deme 1 is never founded by a sampling event from deme 0.  This means that if you simulate for a short period of time, it is unlikely that there is a MRCA common to all individuals sampled.  (See comment above about appropriateness of these examples for research...)_

Usage:

~~~
./migsel_ind N 4Nu_neut 4Nu_sel 4Nr 4Nm s h f1 f2 ngens n outfilename seed
~~~

where:

* N = population number for each deme.
* 4Nu\_neut = 4N times the neutral mutation rate per gamete
* 4Nu\_sel = 4N times the mutation rate per gamete to selected mutations
* 4Nr = 4N times the recombination rate per diploid per region per generation
* 4Nm = 4N times the probability of migration per diploid per generation
* s = the selection coefficient of newly-arising selected mutations
* h = dominance of newly-arising selected mutations
* f1 and f2 are the probabilities of inbreeding in deme 1 and 2, respectively
* ngens = # of generations to simulate
* n = sample size to take from each deme @ end of simulation (must be even #)
* outfilename = name of output file
* seed = random number seed.

Notes: s is taken to be s in deme 1 and -s in deme 2. This lets me illustrate how different fitness functions can be passed to different demes using fwdpp.

The output file contains the following:

1. The metapopulation in binary format
2. Two "ms"-format blocks of size 2*n. These are for neutral and selected mutations, respectively. Within each block, the first n haplotypes are from deme 1 and the second n are from deme 2.

The program writes the data to the output file and then reads it in again. This is mainly to illustrate the binary I/O routines for individual-based metapopulation simulations.

####migsel\_split\_ind (migsel_split_ind.cc) and migsel\_split\_ind\_list (migsel_split_ind_list.cc)

__NOTE:__ migsel_split_ind is not longer compiled as of fwdpp 0.3.0.  It may be resurrected in the future.

Simulates a constant-sized population of N diploids with selection at strength s and dominance h for ngens generations. Then, a copy is made of the population, and the two demes are simulated with migration for another ngens2 generations.

These two programs also serialize their data, read the copied data back, and make sure that input equals output.

These two examples are hints on how to code up "IM"-like models.

Usage:

~~~
./migsel_split_ind N 4Nu_neut 4Nu_sel 4Nr 4Nm s h f1 f2 ngens ngens2 n seed
~~~

~~~
./migsel_split_ind_list N 4Nu_neut 4Nu_sel 4Nr 4Nm s h f1 f2 ngens ngens2 n seed
~~~

where:

* N = population number for each deme.
* 4Nu\_neut = 4N times the neutral mutation rate per gamete
* 4Nu\_sel = 4N times the mutation rate per gamete to selected mutations
* 4Nr = 4N times the recombination rate per diploid per region per generation
* 4Nm = 4N times the probability of migration per diploid per generation
* s = the selection coefficient of newly-arising selected mutations
* h = dominance of newly-arising selected mutations
* f1 and f2 are the probabilities of inbreeding in deme 1 and 2, respectively
* ngens = # of generations to simulate
* ngens2 = # of generations to simulate after the split
* n = sample size to take from each deme @ end of simulation (must be even #)
* seed = random number seed.

Notes: s is taken to be s in deme 1 and -s in deme 2. This lets me illustrate how different fitness functions can be passed to different demes using fwdpp.



####bneck\_selection\_ind (bneck_selection_ind.cc)

This program simulates a population for g generations at size N. In generation g+1, N changes to N2 <= N. The population then grows exponentially to size N3 >= N2 in g2 generations. Selected and neutral mutations are allowed each generation. The output is in “ms” format--one block for neutral mutations followed by one block for selected mutations.

Usage:

~~~
./bneck_selection_ind N theta_neutral theta_sel rho s h g1 N2 N3 g2 n nreps seed
~~~

Where:

* N = starting population size
* theta\_neutral = 4N*(neutral mutation rate per gamete)
* theta\_sel = 4N*(mutation rate per gamete to selected mutations)
* rho = 4Nr, where r is recombination rate per diploid per generation
* s = selection coefficient. Can be negative or positive.
* h = dominance of selected mutations
* N2 = size of bottlenecked population
* N3 = size of recovered population
* g2 = generations taken to go from size N2 to size N3
* n = sample size to take from the population.
* nreps = # replicates to simulate
* seed = random number seed.

####diploid\_ind\_2locus (diploid_ind_2locus.cc)

This program simulates a 2-locus neutral model.  Mutations arise at rate \f$\theta\f$ at each locus, and each locus recombines at rate \f$\rho\f$.  The recombination rate between loci is \f$r_{bw}\f$, which corresponds to the probability a crossover is observed between loci.  In other words, the genotypes at the two loci will switch from \f[\frac{AB}{ab}\f] to \f[\frac{Ab}{aB}\f] with probability \f$r_{bw}\f$.

The output of the simulation is an "ms"-style block for each locus. Positions  are uniform on the interval \f$[0,1)\f$ at locus 1, and uniform on the interval \f$[1,2)\f$ at locus 2.

The usage is:

~~~
Usage:
./diploid_ind_2locus N theta rho rbw ngens n nreps seed
Where:

N = population size (number of diploids)
theta = 4Nu, the scaled neutral mutation rate
rho = 4Nr, the scale recombination rate
rbw = the probability that the two loci cross over, per generation
ngens = the number of generations to simulate
n = the sample size to pull from the population at the end of each simulated replicate
nreps = the number of replicates to simulated
seed = seed value for random number generations
~~~

#### HOC_ind (HOC_ind.cc)

This program simulates a little twist on Turelli's House-of-Cards (HOC) model.  The main point of this model is to demonstrate how to use KTfwd::tags::gamete_dependent to implement mutation models that depend on the current gamete state.

In a typical HOC model (Kingman), the effect size of a mutated allele is a new draw from a Gaussian distribution.  In this simulation, the current effect size of an allele, before mutation, is additive over causative mutations, with sum \f$y\f$.  After the addition of a single mutation, its new effect size will be \f$x\f$, and thus the effect size of the new mutation, \f$e\f$, is chosen such that \f$x = y+e\f$.

The fitness model here is Gaussian stabilizing selection with a mean of 0 and a standard deviation of 1.  With this fitness function and the additive model of Turelli, \f$V(G) \approx 4\mu\f$ when the ratio of the variance in effect sizes over the variance in fitness falls within a certain range.  That approximation seems to hold, at least roughly, for this odd method of generating the effect sizes of new variants.  For example, these commands return a mean of approximately \f$4 \times 10^{-4}\f$

~~~
#Generate 100 seeds
Rscript -e "cat(as.integer(runif(1e2,0,1e6)),sep=\"\n\")" > seeds
#run jobs in parallel on a 4-core machine, write VG to a file
parallel --jobs 4 ./HOC_ind 10000 0.0001 0.25 0 100000 1 {} :::: seeds > VG.txt
~~~

The usage is:

~~~
Usage: ./HOC_ind N mu sigmu rho ngens nreps seed
where:
N = diploid population size
mu = mutation rate to variants affecting fitness
sigmu = standard deviation of effect sizes.  E ~ N(0,sigmu^2)
4Nr = population scaled recombination rate (per diploid, per generation)
ngens = number of generations to simulate
nreps = number of replicates to simulate
seed = random number seed
~~~

The only output is \f$V(G)\f$, the genetic variance in fitness, which is written to stdout for each replicate.

