#!/bin/bash
for i in {1..50} 
do
	sbatch /projects/kernlab/cbattey2/spaceness/batch_scripts/runslim_spatial.srun
done

