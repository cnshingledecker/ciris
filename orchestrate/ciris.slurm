#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:10:00
#SBATCH -p serial
#SBATCH -A herbstgroup

# The program should be in a directory named after the ticket, with the
# params file suitably reflecting the params in the ticket

# START

# 1) Run the program
./ciris

# 2) After the program has run, get the ticket name (current directory)
dirname=(`basename $PWD`)

# 3) Copy params.dat -> ../../done/$dirname
cp params.dat ../../done/$dirname

# 4) Copy slurm output to accounting
cp slurm-* ../../accounting/$dirname

# 5) Remove the current directory
rm -rf ../$dirname

# END
