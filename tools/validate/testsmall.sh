#!bash

parallel --jobs 4 "../../examples/diploid_ind 1000 50 50 10000 10 250 {} > testout.small.{}.txt" ::: $RANDOM $RANDOM $RANDOM $RANDOM

cat testout.small.*.txt | msstats > testout.small.stats

ms 10 1000 -t 50 -r 50 1000 | msstats > testout.small.ms.stats

Rscript plot.R testout.small.stats testout.small.ms.stats testout.small.pdf
