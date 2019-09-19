#!/bin/bash

USAGE="
Usage:
   $0 FILE
where FILE is the name of a SLiM script.
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 0
fi

FILE="$1"

TAG=$(printf "%06d" $RANDOM); 
OUTDIR=$(dirname $FILE)
SLIMFILE=$(basename $FILE)
echo "Directory: $OUTDIR"

LOGFILE=${FILE%.slim}.slurm.log

export OUTDIR
export SLIMFILE
sbatch -o $LOGFILE -e $LOGFILE ./run_sim.sbatch 

