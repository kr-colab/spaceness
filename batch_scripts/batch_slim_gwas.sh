#!/bin/bash
cd ~/projects/spaceness/gwas/sims/
files=./*
for f in $files
do
 
command="Rscript ~/projects/spaceness/demography/smcpp/run_smcpp.R \
--infile $f \
--outdir ~/projects/spaceness/demography/smcpp/output \
--cores 28 \
--recomb 1e-8 \
--mu 1e-8"

cp ~/projects/spaceness/demography/smcpp/smcpp_slurm_header.srun ~/projects/spaceness/demography/smcpp/tmp.srun

echo $command >> ~/projects/spaceness/demography/smcpp/tmp.srun

sbatch ~/projects/spaceness/demography/smcpp/tmp.srun
	
sleep 2

done

