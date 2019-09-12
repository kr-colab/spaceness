#!/bin/bash
for i in {1..200} 
do
	sbatch /projects/kernlab/cbattey2/spaceness/batch_scripts/runslim_random_mating.srun
done

