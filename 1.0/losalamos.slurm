#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 8:00:00
#SBATCH -p serial
#SBATCH -A herbstgroup
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cns7ae@virginia.edu

# Run program
./blastoff
