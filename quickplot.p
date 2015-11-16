#This is a gnuplot script for quickly plotting data from the terminal. 
#It is extremely useful for seeing data right away when you just want
#to do a sanity check on your results.
 set terminal dumb
set logscale x
set datafile separator ','
#Plot vs. fluence
plot 'abundance.csv' using 2:4
#Plot vs. proton count
#plot 'abundance.csv' using 5:4
