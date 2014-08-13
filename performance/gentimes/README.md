When working on fwdpp version 0.2.5, I made a modification of the example program diploid_fixed_sh_ind that timed how long each generation took. (I did this using the C++11 "chrono" header)/

The command line was the following, and I varied run times from 10k to 30k generations, and printed out how long each generation took:
```
./diploid_fixed_sh_ind 1000 400 100 500 -0.001 0. 30000 50 1 101
```

The R script in this folder results in this plot, which shows that run time per generation doesn't really hit an equilibrium until about 10N generations:

![Generation times](gentimes.pdf?raw=true "Time per generation")