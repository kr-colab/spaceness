#!/bin/bash
# get summary statistics for all tree sequences in a directory
cd ~/projects/spaceness

#spatial mate choice sims
files=~/projects/spaceness/sims/msp/*
for f in $files
do
 
 command="python scripts/msp_ts_to_sumstats.py \
 --infile $f \
 --outfile ~/projects/spaceness/sumstats/ss_msp.txt"

cp ~/projects/spaceness/batch_scripts/slurm_sumstats_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.05

done