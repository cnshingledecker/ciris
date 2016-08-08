#!/bin/bash

# update the species.dat file
PARAMS=`cat params.dat`
for line in $PARAMS; do
    if [[ $line =~ ^O,.* ]]; then
        sed -i "s/^O,.*/$line/" species.dat
#    elif [[ $line =~ ^O3,.* ]]; then
#        sed -i "s/^O3,.*/$line/" species.dat
    fi
done
