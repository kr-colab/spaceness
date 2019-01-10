#!/bin/bash
for i in {1..500} 
do
	sbatch /projects/kernlab/cbattey2/spaceness/batch_scripts/runslim_spatial.srun  
	sbatch /projects/kernlab/cbattey2/spaceness/batch_scripts/runslim_random_mating.srun
done

