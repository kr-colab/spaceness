#!/bin/bash
# get summary statistics for all tree sequences in a directory
cd ~/projects/spaceness

#spatial mate choice sims
files=~/projects/spaceness/sims/slimout/spatial/W35/*
for f in $files
do
 
 command="python scripts/sumstats_from_treeseq.py \
 --infile $f \
 --outfile ~/projects/spaceness/sumstats/ss_spatial_W35.txt \
 --sampling random \
 --mu 0.25e-8 \
 --sampling_location 25,25 \
 --seed 12345"

cp ~/projects/spaceness/batch_scripts/slurm_sumstats_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.05

done

#random mating sims
files=~/projects/spaceness/sims/slimout/random_mating/W35/*
for f in $files
do
 
 command="python scripts/sumstats_from_treeseq.py \
 --infile $f \
 --outfile ~/projects/spaceness/sumstats/ss_randmates_W35.txt \
 --sampling random \
 --mu 0.25e-8 \
 --sampling_location 25,25 \
 --seed 12345"

cp ~/projects/spaceness/batch_scripts/slurm_sumstats_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.05

done

