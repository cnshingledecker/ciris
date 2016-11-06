module Rivanna 

using Base.Random.uuid4

include("parameter_writer.jl") # ParameterIO module

# environment parameters
QUEUE = ""
ROOT = ""
ISLANDS = []

# set environment parameters
function initialize(queue, root, islands)
    global QUEUE = queue
    global ROOT = root
    global ISLANDS = islands
    return nothing
end

# Create initial population
function initPop(island, size)
    names = []
    for i = 1:size
        # TODO: replace with ParameterIO.writeSolution
        s = ParameterIO.generateSolution()
        push!(names, ParameterIO.writeSolution(s))
    end
    println("created files")
    for n in names
        submitInput("cp local:$n $ROOT/$island/todo")
        run(`rm $n`)
    end
    println("copied files")
end

# Submit a job
function submitJob(island)
    # pick a file from todo
    file = split(submitInput("ls $ROOT/$island/todo"), "\n")[2]
    # submit that file
    submitInput(string("mv $ROOT/$island/todo/", file, " $ROOT/$island/prog/", file))
    # return the ticket id
    jsdl = generateJSDL(file, island)
    ticket = submitInput(string("qsub $QUEUE local:$jsdl"))[28:63]
    run(`rm $jsdl`)
    return ticket, "$island/prog/$file"
end

# Submits text input to grid command
function submitInput(input)
    if input[end] != '\n'
        input = input * "\n"
    end
    return input
end

# Cats a file given a path relative to root
# NOTE: does NOT split!
function fetchCat(filename)
    return strip(submitInput("cat $ROOT/$filename"))
end

# Gets all done files from an island
function fetchDone(island)
    # grid output has an extra entry - the directory being ls'ed
    return split(strip(submitInput("ls $ROOT/$island/done")), "\n")[2:end]
end

# Gets done size for an island
function fetchDoneSize(island)
    # grid output has an extra entry - the directory being ls'ed
    return length(fetchDone(island)) - 1
end

# Gets todo size for an island
function fetchTodoSize(island)
    # grid output has an extra entry - the directory being ls'ed
    return length(split(submitInput("ls $ROOT/$island/todo"), "\n")) - 1
end

# Deletes a file given a path relative to ROOT
function rmFile(filename)
    submitInput("rm $ROOT/$filename")
end

# Cleans done and prog subdirectories (useful for testing)
function purgeTodo()
    for island in ISLANDS
        submitInput("rm $ROOT/$island/todo/*")
    end
end

function purgeDone()
    for island in ISLANDS
        submitInput("rm $ROOT/$island/done/*")
    end
end
function purgeProg()
    for island in ISLANDS
        submitInput("rm $ROOT/$island/prog/*")
    end
end

# add more candidates, returns number of new files added
function bolster(island)
end

# drop lowest 20% of candidates, returns number of files deleted
function cull(island)
    # get candidates for removal
    files = fetchDone(island)

    # build an array of fitnesses and map to connect score to file
    scores = []
    fitness_map = Dict()
    for filename in files
        try
          file = split(fetchCat("$island/done/$filename"), "\n")
          # get line of fitness score -- should be exactly 1
          fitness = float(filter(line -> contains(line, "FITNESS,"), file)[1][9:end])
          push!(scores, fitness)
          fitness_map[fitness] = filename
        catch error
          if isa(error, BoundsError)
            println("$filename is messed up: removing")
            submitInput("rm -rf $filename")
            println("Bad file removed...")
          end
        end
    end
    sort!(scores)

    # the first 80% are fit
    fit = round(Int, 0.8 * length(scores))
    # the last 20% should be removed
    for score in drop(scores, fit)
        # delete each file
        file = fitness_map[score]
        println("deleting file: $file")
        rmFile("$island/done/$file")
    end
    return length(scores) - fit
end

# creates directory hierarchy
function setupArchipelago()
    if contains(submitInput("cd $ROOT"), "No such file or directory")
        println("ERROR: bad path")
        exit(-1)
    else
        subdirs = ["done" "prog" "todo" "accounting"]
        submitInput("mkdir -p $(join(["$a/$b" for a=ISLANDS, b=subdirs], ' '))")
        submitInput("cd")
    end
end

end
