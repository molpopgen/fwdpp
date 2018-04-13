# Example programs

The best documentation of how to use the library are the example simulations in the examples/ subdirectory of the source code repo.

### Running the examples:

This following list of parameters and their definitions is in common to all of the example programs:

* N = the number of diploid individuals to simulate
* theta = 4Nu, where u is the mutation rate per gamete per generation. For example, theta = 100 means that on average, 2Nu = 50 new mutations are entering the population each generation.
* rho = 4Nr, where r is the recombination rate per diploid per generation
* g = the number of generations to simulate. Often, this should be >= 8N at a minimum
* n = the sample size to draw at the end of the simulation. To match the typical modeling assumptions of population genetics, you should have n << 2N.
* nreps = the number of replicates to simulate
* seed = a random number seed. I will use $RANDOM as a seed, referring to the bash shell method to return a random integer.

Note: familiarity with Hudson's "[ms](http://home.uchicago.edu/~rhudson1)" program is helpful for some of what comes below.

#### diploid_ind (diploid_ind.cc)

The first program, diploid, simulates a Wright-Fisher population with mutation, recombination, and drift. To run it:

~~~
./diploid_ind N theta rho g n nreps seed | gzip > outfile.gz
~~~

The data in outfile.gz will be in the same format as Dick Hudson's "ms" coalescent simulator.

Example:

~~~
./diploid_ind 10000 10 10 100000 50 1 $RANDOM | gzip > test_diploid.gz
~~~


#### diploid\_fixed\_sh\_ind (diploid_fixed_sh\_ind.cc)

This program is similar to diploid_ind, but adds an additional mutation rate (theta\_selected = 4Nu\_s, where u\_s is the mutation rate per gamete per generation to mutations with selection coefficient s) to mutations with selection coefficient s and dominance h. Fitness across sites is multiplicative. The output is in "ms" format--one block for neutral mutations followed by one block for selected mutations.

Usage:

~~~
./dipoid_fixed_sh_ind N theta theta_selected rho s h g n nreps seed | gzip > outfile.gz
~~~

For this program, s can be positive or negative, as can h.

