#!/bin/bash
N=75  #n processes to run in background

cd ~/spaceness

###########spatial sims############
files=~/spaceness/sims/slimout/spatial/W35/*
#random sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_spatial_random_W35.txt \
	 --sampling random \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait

#midpoint sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_spatial_midpoint_W35.txt \
	 --sampling midpoint \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait

#point sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_spatial_point_W35.txt \
	 --sampling point \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait

############### random mating sims ###############
files=~/spaceness/sims/slimout/random_mating/W35/*
#random sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_randmates_random_W35.txt \
	 --sampling random \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait 

#midpoint sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_randmates_midpoint_W35.txt \
	 --sampling midpoint \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait

#point sampling
getsumstats(){
	local f=$1
	python scripts/sumstats_from_treeseq.py \
	 --infile $f \
	 --outfile ~/spaceness/sumstats/ss_randmates_point_W35.txt \
	 --sampling point \
	 --mu 0.25e-8 \
	 --sampling_locs 17,17 \
	 --seed 12345
	sleep 0.05
}
(
for f in $files
	do 
	((i=i%N)); ((i++==0)) && wait
	getsumstats "$f" & 
done
)
i=1
wait


