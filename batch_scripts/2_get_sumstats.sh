#!/bin/bash
# get summary statistics for all tree sequences in a directory
cd ~/projects/spaceness

#spatial mate choice sims
files=~/projects/spaceness/sims/slimout/spatial/W50_run3/*
for f in $files
do
 
 command="python scripts/sumstats_from_treeseq.py \
 --infile $f \
 --outfile ~/projects/spaceness/sumstats/ss_spatial_point_W50.txt \
 --sampling point \
 --mu 1e-8 \
 --sampling_locs 12.5,12.5 \
 --gentimes ~/projects/spaceness/W50sp_gentimes.txt \
 --seed 12345"

cp ~/projects/spaceness/batch_scripts/slurm_sumstats_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.05

done

#random mating sims
files=~/projects/spaceness/sims/slimout/random_mating/W50_run3/*
for f in $files
do
 
 command="python scripts/sumstats_from_treeseq.py \
 --infile $f \
 --outfile ~/projects/spaceness/sumstats/ss_randmates_point_W50.txt \
 --sampling point \
 --mu 1e-8 \
 --sampling_locs 12.5,12.5 \
 --gentimes ~/projects/spaceness/W50rm_gentimes.txt \
 --seed 12345"

cp ~/projects/spaceness/batch_scripts/slurm_sumstats_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.05

done

