module Maintainence

include("parameter_writer.jl")

# is todo too small?
function bolster()
    #println("$ID thinks todo is too small")
    if rand(1:MUTATE_CHANCE) == 1
        #mutate
        println("mutate")
        #file = rand(split(readall(`ls $DONE`), '\n'))
        #new_file = ParameterIO.mutateFile("$DONE/$file")
        #println("made a new file: $new_file")
        #run(`mv $new_file $TODO`)
    else
        #breed
        println("breed")
    end
end

# is done too big?
function cull()
    # drop lowest 20% of performers
    println("remove not implemented yet, need fitness results")

    if rand(1:EXILE_CHANCE) == 1
        # exile
        println("exile not implemented yet")
    elseif rand(1:MUTATE_CHANCE) == 1
        # mutate
        println("mutate")
        #file = rand(split(readall(`ls $DONE`), '\n'))
        #new_file = ParameterIO.mutateFile("$DONE/$file")
        #println("made a new file: $new_file")
        #run(`mv $new_file $TODO`)
        #run(`rm $DONE/$file`)
    else
        # breed
        println("breed")
    end
end

end
