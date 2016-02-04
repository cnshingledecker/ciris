spnum=7
mcount=0
dcount=0
while IFS="," read a b c d e  
do
	if (( $a == $spnum || $b == $spnum )); then
    # Either a or b is the species of interest
		(( dcount=dcount+1 ))
        if (( $a == $spnum && $b == $spnum )); then
        # Both reactants are the one of interest
            (( dcount=dcount+1 ))
        fi 
#	        print "$a+$b=$c+$d+$e" | sed -e 's/^[ \t ]*//' 
	elif (( $c == $spnum || $d == $spnum || $e == $spnum )); then
    # One of the products is
		(( mcount=mcount+1 ))
        if (( $c == $spnum && $d == $spnum && $e == $spnum )); then
            (( mcount=mcount+2 ))
        elif (( $c != $spnum && $d == $spnum && $e == $spnum )); then
            (( mcount=mcount+1 ))
        elif (( $c == $spnum && $d == $spnum && $e != $spnum )); then
            (( mcount=mcount+1 ))
        elif (( $c == $spnum && $d != $spnum && $e == $spnum )); then
            (( mcount=mcount+1 ))
        fi
#	        print "$a+$b=$c+$d+$e" | sed -e 's/^[ \t ]*//'
	fi
done < "reaction_analytics.csv"
print "Species number is $spnum" 
print "The species was made $mcount times and destroyed $dcount times"
print "The difference between creation and destruction is $((mcount-dcount))"
print "The creation-to-destruction ratio is $(( (mcount*1.0)/(dcount*1.0) ))"
