# Example programs

The best documentation of how to use the library are the example simulations in the examples/ subdirectory.

Example programs with _ind in their name are individual-based.  Otherwise, they are gamete-based.

The individual-based programs represent the best/most efficient way to use __fwdpp__.  The others are here for legacy purposes.

###Compiling the examples:

As of fwdpp 0.2.5, the examples are compiled by saying "make check" after the configure step.

###Running the examples:

This following list of parameters and their definitions is in common to all of the example programs:<br>
N = the number of diploid individuals to simulate<br>
theta = 4Nu, where u is the mutation rate per gamete per generation. For example, theta = 100 means that on average, 2Nu = 50 new mutations are entering the population each generation.<br>
rho = 4Nr, where r is the recombination rate per diploid per generation<br>
g = the number of generations to simulate. Often, this should be >= 8N at a minimum<br>
n = the sample size to draw at the end of the simulation. To match the typical modeling assumptions of population genetics, you should have n << 2N.<br>
nreps = the number of replicates to simulate<br>
seed = a random number seed. I will use $RANDOM as a seed, referring to the bash shell method to return a random integer.<br>

Note: familiarity with Hudson's "[ms](http://home.uchicago.edu/~rhudson1)" program is helpful for some of what comes below.

####diploid (diploid.cc)

The first program, diploid, simulates a Wright-Fisher population with mutation, recombination, and drift. To run it:

~~~
./diploid N theta rho g n nreps seed | gzip > outfile.gz
~~~

The data in outfile.gz will be in the same format as Dick Hudson's "ms" coalescent simulator.

Example:

~~~
./diploid 10000 10 10 100000 50 1 $RANDOM | gzip > test_diploid.gz
~~~

####diploid\_ind (diploid_ind.cc)

Identical to diploids, but uses the individual-based sampler. As this program (and diploid) simulate neutral models, this program will typically be slower than that gamete-based implementation when N is large.

####diploid\_binaryIO (diploid_binaryIO.cc)

The next program is called diploid\_binaryIO. This program is identical to diploid, except that it only simulates one replicate at a time, creates two output files, and outputs the entire population rather than a sample of size n << 2N. The first file is an index file, containing an integer representing the replicate number of the output and the position in the haplotypes file where this record begins. The second file is the haplotypes file, which contains the entire population in binary format. 

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
diploid_binaryIO 10000 10 10 100000 $SGE_TASK_ID indexfile hapfile $seed
~~~

The above script, when submitted to a Grid Engine queue, will result in 100 populations of size N=10,000 being written to hapfile. Further, “indexfile” will contain the ID number and position of each file. Records are not over-written because the program uses POSIX file locking to ensure that only 1 process at a time can do the writing. This is a complex program, as it mixes C++ objects with output streams such that they can be written to C-style file descriptors, which is required in order to use file locking (which is a C feature with no C++ analog). However, the advantage is that you write all data to one large file, avoiding the plague of lots of small files that can bring distributed file systems to their knees.

####diploid\_binaryIO\_ind (diploid_binaryIO_ind.cc)

Identical to diploid\_binaryIO, but individual-based.

####diploid\_fixed\_sh (diploid_fixed_sh.cc)

This program is similar to diploid, but adds an additional mutation rate (theta\_selected = 4Nu\_s, where u\_s is the mutation rate per gamete per generation to mutations with selection coefficient s) to mutations with selection coefficient s and dominance h. Fitness across sites is multiplicative. The output is in "ms" format--one block for neutral mutations followed by one block for selected mutations.

Usage:

~~~
./dipoid_fixed_sh N theta theta_selected rho s h g n nreps seed | gzip > outfile.gz
~~~

For this program, s can be positive or negative, as can h.

####diploid\_fixed\_sh\_ind (diploid_fixed_sh_ind.cc) and diploid\_fixed\_sh\_ind_lambda (diploid_fixed_sh_ind_lambda.cc)

Identical to diploid\_fixed\_sh, but based on the individual-based sampler.

These two programs are identical in function, but differ in that the former is implemented using std::bind and the latter using lambda expressions.

####diploid\_twopop\_mig (diploid_twopop_mig.cc)

This program simulates an ancestral population for g generations, at which point a daughter population of size N is “budded” off from the ancestral population. Evolution continues for g2 more generations, with symmetric migration at rate M = 4Nm, where m is the migration rate per diploid per generation. The output is in "ms" format, with n haplotypes per sample just like how ms outputs data for multi-population models.

Usage:

~~~
./diploid_twopop_mig N theta rho g g2 M n nreps seed | gzip > outfile.gz
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
N = population number for each deme.
4Nu\_neut = 4N times the neutral mutation rate per gamete
4Nu\_sel = 4N times the mutation rate per gamete to selected mutations
4Nr = 4N times the recombination rate per diploid per region per generation
4Nm = 4N times the probability of migration per diploid per generation
s = the selection coefficient of newly-arising selected mutations
h = dominance of newly-arising selected mutations
f1 and f2 are the probabilities of inbreeding in deme 1 and 2, respectively
ngens = # of generations to simulate
n = sample size to take from each deme @ end of simulation (must be even #)
outfilename = name of output file
seed = random number seed.

Notes: s is taken to be s in deme 1 and -s in deme 2. This lets me illustrate how different fitness functions can be passed to different demes using fwdpp.

The output file contains the following:

1. The metapopulation in binary format
2. Two "ms"-format blocks of size 2*n. These are for neutral and selected mutations, respectively. Within each block, the first n haplotypes are from deme 1 and the second n are from deme 2.

The program writes the data to the output file and then reads it in again. This is mainly to illustrate the binary I/O routines for individual-based metapopulation simulations.

####migsel\_split\_ind (migsel_split_ind.cc) and migsel\_split\_ind\_list (migsel_split_ind_list.cc)

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
N = population number for each deme.
4Nu\_neut = 4N times the neutral mutation rate per gamete
4Nu\_sel = 4N times the mutation rate per gamete to selected mutations
4Nr = 4N times the recombination rate per diploid per region per generation
4Nm = 4N times the probability of migration per diploid per generation
s = the selection coefficient of newly-arising selected mutations
h = dominance of newly-arising selected mutations
f1 and f2 are the probabilities of inbreeding in deme 1 and 2, respectively
ngens = # of generations to simulate
ngens2 = # of generations to simulate after the split
n = sample size to take from each deme @ end of simulation (must be even #)
seed = random number seed.

Notes: s is taken to be s in deme 1 and -s in deme 2. This lets me illustrate how different fitness functions can be passed to different demes using fwdpp.

####RHH (RHH.cc)

This program simulates the "recurrent hitch-hiking" model in which the expected number of fixations of beneficial codominant mutations per site per 4N generations is Lambda. The output is ms-format (just one block for neutral mutations). This model is one of the most thoroughly-studied models of the effect of selection on linked variation. See, for example:

1. Wiehe, T. H. and W. Stephan, 1993 Analysis of a genetic hitchhiking model, and its application to DNA polymorphism data from Drosophila melanogaster. Molecular Biology and Evolution 10: 842–854.
2. Kaplan, N. L., R. R. Hudson and C. H. Langley, 1989 The "hitchhiking effect" revisited. Genetics 123: 887–899.
4. Braverman, J. M., R. R. Hudson, N. L. Kaplan, C. H. Langley and W. Stephan, 1995 The hitchhiking effect on the site frequency spectrum of DNA polymorphisms. Genetics 140: 783–796.
4. Przeworski, M., 2002 The Signature of Positive Selection at Randomly Chosen Loci. Genetics 160: 1179.

The model in this program allows selected mutations within the region with neutral mutation rate theta (the so-called "neutral" or "sampled" region in the coalescent literature treating RHH via the structured coalescent). In general, the majority of selected sites will appear up to s/r\_bp base pairs, where r is the recombination rate per base pair.

The population is simulated for g generations with only neutral mutations, followed by g2 more generations with neutral and selected mutations. The selection coefficient, s, must be > 0, but the program doesn't check for this, so just enter an appropriate value please. The output is in ms format, and only contains data for neutral mutations, so any segregating selected mutations will not be output. The region with mutation rate theta is L base pairs long. L is really only used to calculate r\_bp, which is rho/(4*N*(L-1)), as there are L-1 possible positions for recombination to occur in a region of L nucleotides.

Usage:

~~~
./RHH N theta rho L s Lambda g1 g2 n nreps seed | gzip > outfile.gz
~~~

Note: many of the well-known formulas for the effect of RHH on linked, neutral variation make very strong assumptions about the parameter values. For example, N needs to be large and s needs to be small, but Ns needs to be large. Further, Lambda needs to be sufficiently small such that sweeps are independent in time. This means that plugging in values to this program and comparing to theoretical predictions may lead to apparent discrepancies. This is also the case with various coalescent simulations of RHH.

####bneck\_selection\_ind (bneck_selection_ind.cc)

This program simulates a population for g generations at size N. In generation g+1, N changes to N2 <= N. The population then grows exponentially to size N3 >= N2 in g2 generations. Selected and neutral mutations are allowed each generation. The output is in “ms” format--one block for neutral mutations followed by one block for selected mutations.

Usage:

~~~
./bneck_selection_ind N theta_neutral theta_sel rho s h g1 N2 N3 g2 n nreps seed
~~~

Where:<br>
N = starting population size<br>
theta\_neutral = 4N*(neutral mutation rate per gamete)<br>
theta\_sel = 4N*(mutation rate per gamete to selected mutations)<br>
rho = 4Nr, where r is recombination rate per diploid per generation<br>
s = selection coefficient. Can be negative or positive.<br>
h = dominance of selected mutations<br>
N2 = size of bottlenecked population<br>
N3 = size of recovered population<br>
g2 = generations taken to go from size N2 to size N3<br>
n = sample size to take from the population.<br>
nreps = # replicates to simulate<br>
seed = random number seed.<br>

####TFL2013 (TFL2013.cc)

This is the simulation from Thornton, Foran, and Long (2013) Properties and Modeling of GWAS when Complex Disease Risk Is Due to Non-Complementing, Deleterious Mutations in Genes of Large Effect. PLoS Genetics 9: e1003258. The manuscript is available [here](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003258). See that paper for modeling details.

Usage:

~~~
./TFL2013 N mu_disease mu_neutral littler esize dist sdE sdS burnin evolve msfile phenofile effectfile
~~~

Where:<br>
N = diploid population number<br>
mu\_disease = mutation rate to causative mutations (per gamete)<br>
mu\_neutral = mutation rate to neutral mutations (per gamete)<br>
littler = Recombination rate per diploid<br>
esize = Effect size of causative mutation.<br>
dist = boolean.  If 0, then esize is fixed.  If 1, then esize is the mean of an exponentural with mean esize<br>
sdE = standard deviation in Gaussian noise added to phenotype<br>
sdS = standard deviation in Gaussian fitness function. (Should probably just set to 1.)<br>
burnin = number of generations to evolve with no causative mutations.<br>
evolve = number of generations to evolve with neutral and causative mutations.<br>
msfile = File to write the entire population in "ms" format.  Each sequential pair of gametes is a different diploids. Output is gzipped.<br>
phenofile = File to write phenotypes of population.  Two columns: genetic component first, then Gaussian noise component.  The phenotype is the sum of the two.  Order of phenotypes is same as order of diploids in msfile.  Gzipped<br>
effectfile = File to write positions and effect sizes of causative mutations. Gzipped.<br>

Example from the PLoS Genetics paper:

~~~
./TFL2013 20000 0.000125 0.00125 0.00125 0.5 1 0.075 1 0 160000 msfile.gz phenotypes.gz effects.gz
~~~

will simulate a mean effect size of 0.5 (exponentially-distributed).  The above will take probably many hours to run.

####diploid\_ind\_2locus (diploid_ind_2locus.cc)

This program simulates a 2-locus neutral model.  Mutations arise at rate \f$\theta\f$ at each locus, and each locus recombines at rate \f$\rho\f$.  The recombination rate between loci is \f$r_{bw}\f$, which corresponds to the probability a crossover is observed between loci.  In other words, the genotypes at the two loci will switch from \f[\frac{AB}{ab}\f] to \f[\frac{Ab}{aB}\f] with probability \f$r_{bw}\f$.

The output of the simulation is an "ms"-style block for each locus. Positions  are uniform on the interval \f$[0,1)\f$ at locus 1, and uniform on the interval \f$[1,2)\f$ at locus 2.

The usage is:

~~~
Usage:
./diploid_ind_2locus N theta rho ngens n nreps seed
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
r = recombination rate (per diploid, per generation)
nreps = number of replicates to simulate
seed = random number seed
~~~

The only output is \f$V(G)\f$, the genetic variance in fitness, which is written to stdout for each replicate.

