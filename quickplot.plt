#This is a gnuplot script for quickly plotting data from the terminal. 
#It is extremely useful for seeing data right away when you just want
#to do a sanity check on your results.
#set terminal postscript eps enhanced color font '~/hershey-fonts/hershey-fonts/futural.jhf'
#set output 'plotabundance.eps'
#set terminal dumb 
set term x11
set logscale x
#set logscale y
#set datafile separator ' '
#Plot vs. fluence
#set term x11 0
plot 'abundance.wsv' using 1:3 with line, \
     'abundance.wsv' using 1:5 with line
#set term x11 1
#plot 'abundance.csv' using 2:6 with line
pause -1
#Plot vs. proton count
#plot 'abundance.csv' using 5:4

