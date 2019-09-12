#!/bin/bash
for i in {1..200} 
do
	sbatch /projects/kernlab/cbattey2/spaceness/scripts/runslim_wf.srun
done

