#!/bin/bash

# need to grab random seed used
if [[ -e "seed" ]]; then
    seed=`cat seed`
    echo "SEED,$seed" >> "params.dat"
fi

# need to grab fitness score
if [[ -e "fitness_results" ]]; then
    # grabs the last line
    tokens=( `tail -1 fitness_results` )
    # grabs the fitness score - last element of 4-tuple
    fitness=${tokens[3]}
    # $CANDIDATE should be a defined shell variable
    # contains the name of the candidate solution's file
    echo "FITNESS,$fitness" >> "params.dat"
fi
