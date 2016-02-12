#This is a gnuplot script for quickly plotting data from the terminal. 
#It is extremely useful for seeing data right away when you just want
#to do a sanity check on your results.
#set terminal postscript eps enhanced color font '~/hershey-fonts/hershey-fonts/futural.jhf'
#set output 'plotabundance.eps'
#set terminal dumb 
set logscale x
set datafile separator ','
#Plot vs. fluence
plot 'abundance.csv' using 2:4 with line
#'abundance.csv' using 2:5 with line
pause -1
#Plot vs. proton count
#plot 'abundance.csv' using 5:4

