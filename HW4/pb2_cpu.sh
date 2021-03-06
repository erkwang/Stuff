#!/bin/bash -l

module load R/3.0.0

#SBATCH --job-name=pb2_cpu
#SBATCH --output=pb2_cpu.out
#SBATCH --error=pb2_cpu.err

srun R --no-save --vanilla < pb2_cpu.R
