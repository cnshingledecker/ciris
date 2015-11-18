#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 20:00:00
#SBATCH -p serial
#SBATCH -A herbstgroup

# Run Program
module load gcc/5.2.0
python vary_params.py
