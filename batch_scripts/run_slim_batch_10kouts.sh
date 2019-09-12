#!/bin/bash
for i in {1..100} 
do
	sbatch /projects/kernlab/cbattey2/spaceness/scripts/runslim_10kouts.srun
done

