#Using C++11 in fwdpp

Veresion 0.2.5 of fwdpp is the beginning of a migration from "classic" C++ to the C++11 language standard.   In practice, that means that many very useful features that have only been available via the [boost](http://www.boost.org) libraries are now available in the C++ standard library.

This directory represents a set of tests that I performed on the UCI cluster to check that "de-boosting" fwdpp won't adversely affect performance.

##Test1: std::bind vs. boost::bind

This test is of the following conditions:

1. The dev branch of fwdpp (which is the same as the published version 0.2.4 as far as this is concerned) was compiled with -std=c++11 using GCC 4.7.3 on the UCI HPC, which is a cluster of AMD0-based 64-core nodes.  The example programs were compiled against libsequence 1.8.0 and boost 1.54.0
2. The cpp11 branch of fwdpp where the only difference from the above is that boost::bind was replace with std::bind
3. The example program diploid_fixed_sh_ind was used to simulate a population of N=10,000 diploids for 8N generations with 4Nu = 4Nr and 4Nd = 5, which is twice the number of deleterious mutations entering the population per generation).  The script to run the simulation is:

```{sh}
#!/bin/bash

cd $SGE_O_WORKDIR

module load krthornt/libsequence/1.8.0
module load boost/1.54.0

if [ ! -s time.$SGE_TASK_ID.txt ]
then
    /usr/bin/time -f "%e %M" -o time.$SGE_TASK_ID.txt $1/diploid_fixed_sh_ind 10000 $2 5 $3 -0.001 0.5 80000 50 1 $SGE_TASK_ID > out.$SGE_TASK_ID.txt
fi
```

The run times and memory usage are basically the same using std::bind vs boost::bind, with an arguably advantage to std::bind.  In the figures, "0.2.4" refers to the dev branch/published version of fwdpp, and "0.2.5" will refer to the C++11 branch.

Run times:
![test1time](t1_time.png?raw=true "Run time and bind")

Memory (you only see one line because there is no difference between the two versions, so they overlap completely):
![test1mem](t1_mem.png?raw=true "Memory and bind")

##TEST 2: flly-cpp11 + fewer recombination searches + GCC 4.8