#!/bin/bash
for sigma in `seq 0.2 .25 3`
do
echo $sigma
slim -t -d sigma=$sigma -d outpath="'~/Desktop/test'" slim_recipes/spaceness.slim | grep "CPU time used"
done

