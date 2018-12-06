#!/bin/bash
# get summary statistics for all tree sequences in a directory
cd ~/projects/spaceness

#spatial mate choice sims
files=~/projects/spaceness/sims/slimout/spatial/W35/*
for f in $files
do
 
 command="python scripts/run_gwas.py \
 --treeseq $f \
 --outdir ~/projects/spaceness/gwas/out \
 --plink_path ~/projects/plink-1.07-x86_64/plink \
 --vcftools_path ~/projects/vcftools_0.1.13/bin/vcftools \
 --nSamples 1000 \
 --sampling random \
 --mu 0.25e-8 \
 --phenotype transform_coord \
 --phenotype_mean 175 \
 --phenotype_sd 8 \
 --seed 12345"
 
cp ~/projects/spaceness/batch_scripts/slurm_gwas_header.srun ~/projects/spaceness/batch_scripts/tmp.srun

echo $command >> ~/projects/spaceness/batch_scripts/tmp.srun

sbatch ~/projects/spaceness/batch_scripts/tmp.srun
	
sleep 0.01

done
