#!/bin/bash

for i in {1..10}
do 
	simnum=$i
	command="slim -d sigma=0.2 -d simnum=$simnum /projects/kernlab/cbattey2/spaceness/slim_recipes/flat_map_vcfout.slim"
	cp ~/projects/spaceness/batch_scripts/runslim_header.srun ~/projects/spaceness/batch_scripts/tmp.srun
	echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun
	#tail ~/projects/spaceness/batch_scripts/tmp.srun
	sbatch ~/projects/spaceness/batch_scripts/tmp.srun
done