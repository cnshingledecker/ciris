#! /bin/bash

rm elec_nums.csv 
rm abundance.csv
rm test_wait_list.txt
rm wait_list_flaw.txt
echo "Removed old output files"
echo "Starting BLASTOFF"
./blastoff
