#!bash

parallel --jobs 4 "../../examples/diploid_ind 5000 500 500 50000 100 25 {} > testout.large.{}.txt" ::: $RANDOM $RANDOM $RANDOM $RANDOM

cat testout.large.*.txt | msstats > testout.large.stats

ms 100 1000 -t 500 -r 500 1000 | msstats > testout.large.ms.stats

Rscript plot.R testout.large.stats testout.large.ms.stats testout.large.pdf
