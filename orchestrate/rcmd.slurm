#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 24:00:00
#SBATCH -p serial
#SBATCH -A herbstgroup

# The program should be in a directory named after the ticket, with the
# params file suitably reflecting the params in the ticket

module load julia

# START
./RCMD
# END
