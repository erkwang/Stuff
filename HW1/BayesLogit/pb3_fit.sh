#!/bin/bash -l

module load R/3.0.0

#SBATCH --job-name=pb3_fit
#SBATCH --output=dump/pb3_fit.out
#SBATCH --error=dump/pb3_fit.err

srun R --no-save --vanilla < pb3_fit.R



