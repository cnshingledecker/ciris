module Grid

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

# initializes grid process
function initGrid()
    print("Starting grid process...\t")
    g_out, g_in, g_proc = readandwrite(`/home/alex/GenesisII/bin/grid`)
    readuntil(g_out, "[grid] ")
    println("Done.")
    return g_out, g_in, g_proc
end

# references to grid process IO and status
g_out, g_in, g_proc = initGrid()

# read in JSDL file, sub placeholders, and write out JSDL
function generateJSDL(target_file, island)
    jsdl = readall("losalamos.jsdl")
    jsdl = replace(jsdl, "PLACEHOLDER0", target_file)
    jsdl = replace(jsdl, "PLACEHOLDER1", island)
    file = open("tmp.jsdl", "w")
    write(file, jsdl)
    close(file)
    return "tmp.jsdl"
end

# Create initial population
function initPop(island, size)
    names = []
    for i = 1:size
        # TODO: replace with ParameterIO.writeSolution
        push!(names, string(uuid4()))
        writecsv(names[i], ParameterIO.generateSolution())
    end
    local_names = map(n -> "local:$n", names)
    println("created files")
    submitInput("cp $(join(local_names, ' ')) $ROOT/$island/todo")
    println("copied files")
    for n = names
        run(`rm $n`)
    end
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
    if !process_running(Grid.g_proc)
        Grid.g_out, Grid.g_in, Grid.g_proc = initGrid()
    end
    if input[end] != '\n'
        input = input * "\n"
    end
    write(Grid.g_in, input)
    response = strip(readuntil(Grid.g_out, "[grid] "))
    return replace(response, "[grid]", "")
end

# Gets all jobs running in queue
function fetchJobQueue()
    return split(strip(submitInput("qstat $QUEUE")), "\n")[2:end]
end

# Gets a list of all tickets in the queue
function fetchQueueTickets()
    return map(z -> first(split(z)), fetchJobQueue())
end

# Gets a list of all tickets in the queue with status FINISHED or FAILED
function fetchCompletedTickets()
    jobs = filter(line -> contains(line, "FINISHED") || contains(line, "FAILED"), fetchJobQueue())
    return map(z -> first(split(z)), jobs)
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

# Copies a local file to the grid relative to ROOT
function cpLG(local_path, grid_path)
    submitInput("cp local:$local_path $ROOT/$grid_path")
end

# Removes finished job from queue and file given a path relative to ROOT
function rmJob(ticket, filename)
    submitInput("qcomplete $QUEUE $ticket")
    submitInput("rm $ROOT/$filename")
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
        file = split(fetchCat("$island/done/$filename"), "\n")
        # get line of fitness score -- should be exactly 1
        fitness = float(filter(line -> contains(line, "FITNESS,"), file)[1][9:end])
        push!(scores, fitness)
        fitness_map[fitness] = filename
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
    if contains(submitInput("cd $ROOT"), "does not exist")
        println("[Grid.setupArchipelago] ERROR: bad path")
        exit(-1)
    else
        subdirs = ["done" "prog" "todo" "accounting"]
        submitInput("mkdir $(join(ISLANDS, ' '))")
        submitInput("mkdir $(join(["$a/$b" for a=ISLANDS, b=subdirs], ' '))")
        submitInput("cd")
    end
end

function cleanup()
    if !process_running(Grid.g_proc)
        submitInput("exit")
        close(Grid.g_in)
        close(Grid.g_out)
    end
    quit()
end

end
