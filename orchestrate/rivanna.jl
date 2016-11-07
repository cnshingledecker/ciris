#=
For full Rivanna compatability, remove occurances of:
1) local
2) queue

Also check
1) submit job returns the uuid of the job submitted
=#
module Rivanna

using Base.Random.uuid4
include("parameter_writer.jl") # ParameterIO module

# environment parameters
ROOT = ""
ISLANDS = []

# set environment parameters
function initialize(root, islands)
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
        run(`cp $n $ROOT/$island/todo`)
        run(`rm $n`)
    end
    println("copied files")
end

# Submit a job
function submitJob(island)
    # pick first file from todo
    file = readdir("$ROOT/$island/todo") 
    file = file[1]
    # Make a directory in the root directory with the filename
    println("Now making $island/prog/$file")
    run(`mkdir $ROOT/$island/prog/$file`)
    cd("$ROOT/$island/prog/$file")
    println("now in directory",pwd())
    # Populate the new directory with the inputs/src
    println("Now copying inputs to $island/prog/$file")
    srcfiles = readdir("$ROOT/src")
    for srcfile in srcfiles
      run(`cp $ROOT/src/$srcfile $ROOT/$island/prog/$file/`)
    end
    run(`chmod 777 $ROOT/$island/prog/$file/ciris.slurm`)
    # Copy the params from TODO to local params.dat
    run(`mv $ROOT/$island/todo/$file $ROOT/$island/prog/$file/params.dat`)
    # Now the new directory with all the name=$file should have all the
    # necessary input: submit slurm script
#    run(`sbatch $ROOT/$island/prog/$file/ciris.slurm`)
    run(`$ROOT/$island/prog/$file/ciris.slurm`)
    # return the ticket id
    return 
end

# Copies a local file to the grid relative to ROOT
function cpLG(local_path, rivanna_path)
    run(`cp $local_path $ROOT/$rivanna_path`)
end

# Cats a file given a path relative to root
# NOTE: does NOT split!
function fetchCat(filename)
    return strip(run(`cat $ROOT/$filename`))
end

# Gets all done files from an island
function fetchDone(island)
    # grid output has an extra entry - the directory being ls'ed
    return split(strip(run(`ls $ROOT/$island/done`)), "\n")[2:end]
end

# Gets done size for an island
function fetchDoneSize(island)
    # grid output has an extra entry - the directory being ls'ed
    length = round(Int,size(readdir("$ROOT/$island/done"),1))
    println("$island/done has $length elements")
    return length
end

# Gets todo size for an island
function fetchTodoSize(island)
    # grid output has an extra entry - the directory being ls'ed
    length = round(Int,size(readdir("$ROOT/$island/todo"),1))
    println("$island/todo has $length elements")
    return length
end

# Deletes a file given a path relative to ROOT
function rmFile(filename)
    run(`rm $ROOT/$filename`)
end

# Cleans done and prog subdirectories (useful for testing)
function purgeTodo()
    for island in ISLANDS
        run(`rm $ROOT/$island/todo/*`)
    end
end

function purgeDone()
    for island in ISLANDS
        run(`rm $ROOT/$island/done/*`)
    end
end
function purgeProg()
    for island in ISLANDS
        run(`rm $ROOT/$island/prog/*`)
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
            run(`rm -rf $filename`)
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
    println("Now in directory",pwd())
    subdirs = ["done" "prog" "todo" "accounting"]
    for island in ISLANDS
        println("Now making directory for $island")
        run(`mkdir -v $ROOT/$island`)
        for subdir in subdirs
            run(`mkdir -v $ROOT/$island/$subdir`)
        end
    end
    cd("$ROOT")
    run(`ls`)
end

end
