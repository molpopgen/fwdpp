#!bash

parallel --jobs 4 "../../examples/diploid_ind 1000 250 250 10000 10 125 {} > testout.small.{}.txt" ::: $RANDOM $RANDOM $RANDOM $RANDOM

cat testout.small.*.txt | msstats > testout.small.stats

ms 10 500 -t 250 -r 250 1000 | msstats > testout.small.ms.stats

Rscript plot.R testout.small.stats testout.small.ms.stats testout.small.pdf
