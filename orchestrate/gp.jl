#!/bin/julia
#=
 = gp.jl
 =  Orchestration program for genetic programming of losalamos
 =
 =
 = Alex Thomas
 =#

include("parameter_writer.jl")
include("grid.jl")
include("maintainence.jl")

# island population constants
MIN_WORK = 3   # max number of jobs/island
MIN_TODO = 50   # min todo size/island
MAX_DONE = 375  # max done size/island
EXILE = 25
MUTATE_CHANCE=10 # 1 out of...
EXILE_CHANCE=5   # 1 out of...

# Some sanity assertions
@assert MIN_WORK < MIN_TODO "Max possible jobs/island exceeds minimum population threshold!"

# reads a ticket store into memory
function readTickets(filename)
    file = readcsv(filename)
    return Dict(zip(file[:,1], file[:,2]))
end

# writes a ticket store into memory
function writeTickets(filename, tickets)
    if isempty(tickets)
        run(`rm $filename`)
    else
        writecsv(filename, tickets)
    end
end

# counts the number of tickets per island
function countTickets(tickets)
    count = Dict()
    for k in keys(tickets)
        # values are of the form <island>/prog/<filename>
        island = split(tickets[k], '/')[1]
        if !haskey(count, island)
            count[island] = 1
        else
            count[island] = count[island] + 1
        end
    end
    return count
end

# name of ticket store (present from previous runs
TICKET_MASTER = "ticket_master.csv"
tickets = isfile(TICKET_MASTER) ? readTickets(TICKET_MASTER) : Dict()

# control logic for maintaining jobs and islands
function do_maintainence(islands)
    # check if todo population should be bolstered
    for island in islands
        todo_size = Grid.fetchTodoSize(island)
        if todo_size < MIN_TODO
            deficit = MIN_TODO - todo_size
            println("$island should be bolstered")

            # add more candidates to todo from done
            done_files = Grid.fetchDone(island)
            if length(done_files) < 2
                println("$island: not enough files in done to bolster")
                continue
            end
            # add "deficit" number of files
            mutated = 0
            bred = 0
            for j=1:deficit
                if rand(1:MUTATE_CHANCE) == 1
                    #mutate
                    # pick a random file and fetch its contents
                    file_contents = Grid.fetchCat("$island/done/$(rand(done_files))")
                    parsed_file = ParameterIO.readSolution(IOBuffer(file_contents))
                    # drop fitness
                    delete!(parsed_file, "FITNESS")
                    # vary number of mutations
                    for i = rand(1:4)
                        ParameterIO.mutate!(parsed_file)
                    end
                    name = ParameterIO.writeSolution(parsed_file)
                    Grid.cpLG(name, "$island/todo")
                    run(`rm $name`)
                    mutated = mutated + 1
                else
                    #breed
                    # pick two files at random
                    f1, f2 = rand(done_files, 2)
                    f1_contents = Grid.fetchCat("$island/done/$f1")
                    f2_contents = Grid.fetchCat("$island/done/$f2")
                    f1_parsed = ParameterIO.readSolution(IOBuffer(f1_contents))
                    f2_parsed = ParameterIO.readSolution(IOBuffer(f2_contents))
                    delete!(f1_parsed, "FITNESS")
                    delete!(f2_parsed, "FITNESS")
                    # breed and stash new candidate
                    child = ParameterIO.breed(f1_parsed, f2_parsed)
                    name = ParameterIO.writeSolution(child)
                    Grid.cpLG(name, "$island/todo")
                    run(`rm $name`)
                    bred = bred + 1
                end
            end
            println("$island: mutated $mutated, bred $bred")
        end
    end

    # clean the queue, prog, and tickets map
    # TODO: if job manually removed, update tickets map
    completed_jobs = Grid.fetchCompletedTickets()
    for ticket in completed_jobs
        Grid.rmJob(ticket, tickets[ticket])
        delete!(tickets, ticket)
    end

    # check if more jobs should be submitted
    ticket_count = countTickets(tickets)
    println("ticket count: ", ticket_count)
    for island in islands
        jobs_in_progress = get(ticket_count, island, 0)
        println("$island jobs: $jobs_in_progress")
        if jobs_in_progress < MIN_WORK
            # below threshold!
            new_jobs = MIN_WORK - jobs_in_progress
            println("Making $new_jobs for $island")
            for j=1:new_jobs
                ticket, path = Grid.submitJob(island)
                tickets[ticket] = path
            end
        end
    end

    # check if done population is too big
    for island in islands
        done_size = Grid.fetchDoneSize(island)
        if done_size > MAX_DONE
            gap = MAX_DONE - done_size
            deleted = Grid.cull(island)
            println("$island: culled $deleted")
        end
    end
end


println("GP SYSTEM FOR LOSALAMOS\n")

if length(ARGS) == 1 && (ARGS[1] == "-d" || ARGS[1] == "--debug")
    # debug/default
    println("Using default arguments")
    queue = "/resources/xcg.virginia.edu/queues/vm-queue"
    root = "/home/xcg.virginia.edu/apt9jf/archipelago"

    islands = ["G0"]
    Grid.initialize(queue, root, islands)

    println("Queue: $(Grid.QUEUE)")
    println("ROOT: $(Grid.ROOT)")
    println("ISLANDS: $(Grid.ISLANDS)")
elseif length(ARGS) >= 3
    Grid.initialize(ARGS[1], ARGS[2], ARGS[3:end])
    println("Queue: $(Grid.QUEUE)")
    println("ROOT: $(Grid.ROOT)")
    println("ISLANDS: $(Grid.ISLANDS)")
else
    # bad
    println("Bad call.")
    println("<path_to_queue> <path_to_root> island0 island1 island2 ...")
    exit(-1)
end

# the main loop
do_maintainence(Grid.ISLANDS)
#while true
#    for i=1:5
#        do_maintainence(Grid.ISLANDS)
#        # sleep 2 seconds
#        sleep(2)
#    end
#    println("Time to quit?")
#    if !isfile("GO")
#        println("Quitting...")
#        break
#    end
#end
writecsv(TICKET_MASTER, tickets)
Grid.cleanup()
