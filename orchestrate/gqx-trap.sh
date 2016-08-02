function finish {
    # re-start service
   ./gp.jl "/resources/xcg.virginia.edu/queues/grid-queue-xcg3" "/home/xcg.virginia.edu/cns7ae/archipelago" "G0" "G1" "G2" "G3" "G4"
}
while :
    do
        trap finish EXIT
        ./gp.jl "/resources/xcg.virginia.edu/queues/grid-queue-xcg3" "/home/xcg.virginia.edu/cns7ae/archipelago" "G0" "G1" "G2" "G3" "G4"
done

