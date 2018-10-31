#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks-per-core=1

module use /projects/apps/shared/modulefiles
module load SLiM

cd $SLURM_SUBMIT_DIR

: ${OUTDIR?Must define OUTDIR}
if [ -z "$OUTDIR" ]
then
    echo "Must define OUTDIR (is empty)."
    exit 1
fi

: ${SLIMFILE?Must define SLIMFILE}
if [ -z "$SLIMFILE" ]
then
    echo "Must define SLIMFILE (is empty)."
    exit 1
fi

LOGFILE=${SLIMFILE}.log

echo "Running: in $OUTDIR"
echo "   slim -d $SLIMFILE \&\> $LOGFILE"

cd $OUTDIR

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    slim $SLIMFILE &>> $LOGFILE

echo "Done!"


