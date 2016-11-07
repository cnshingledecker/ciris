#!/usr/bin/julia
#=
= rgp.jl
=  Orchestration program for genetic programming of CIRIS
=
= Christopher Shingledecker & Alex Thomas
=#

include("parameter_writer.jl")
include("rivanna.jl")
include("maintainence.jl")
using DataFrames

# island population constants
MIN_WORK      = 1  # max number of jobs/island
MIN_TODO      = 10 # min todo size/island
MAX_DONE      = 5 # max done size/island
EXILE         = 25 # not implemented yet...
MUTATE_CHANCE = 5  # 1 out of...
EXILE_CHANCE  = 5  # 1 out of...

# Some sanity assertions
@assert MIN_WORK < MIN_TODO "Max possible jobs/island exceeds minimum population threshold!"


# control logic for maintaining jobs and islands
function do_maintainence(islands)
    # check if todo population should be bolstered
    for island in islands
        todo_size = Rivanna.fetchTodoSize(island)
        println("In rgp.jl, called todo_size in do_maintainence")
        if todo_size < MIN_TODO
            deficit = MIN_TODO - todo_size
            println("$island should be bolstered")
            # add more candidates to todo from done
            done_files = readdir("$(Rivanna.ROOT)/$island/done")
            println("done_files_length=",size(done_files,1))
            if size(done_files,1) < 2
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
                    file_contents = Rivanna.fetchCat("$island/done/$(rand(done_files))")
                    parsed_file = ParameterIO.readSolution(IOBuffer(file_contents))
                    # drop fitness and seed
                    delete!(parsed_file, "FITNESS")
                    delete!(parsed_file, "SEED")
                    # vary number of mutations
                    for i = rand(1:4)
                        ParameterIO.mutate!(parsed_file)
                    end
                    name = ParameterIO.writeSolution(parsed_file)
                    Rivanna.cpLG(name, "$island/todo")
                    run(`rm $name`)
                    mutated = mutated + 1
                else
                    #breed
                    # pick two files at random
                    f1, f2 = rand(done_files, 2)
                    f1_contents = Rivanna.fetchCat("$island/done/$f1")
                    f2_contents = Rivanna.fetchCat("$island/done/$f2")
                    f1_parsed = ParameterIO.readSolution(IOBuffer(f1_contents))
                    f2_parsed = ParameterIO.readSolution(IOBuffer(f2_contents))
                    delete!(f1_parsed, "FITNESS")
                    delete!(f2_parsed, "FITNESS")
                    delete!(f1_parsed, "SEED")
                    delete!(f2_parsed, "SEED")
                    # breed and stash new candidate
                    child = ParameterIO.breed(f1_parsed, f2_parsed)
                    name = ParameterIO.writeSolution(child)
                    Rivanna.cpLG(name, "$island/todo")
                    run(`rm $name`)
                    bred = bred + 1
                end
            end
            println("$island: mutated $mutated, bred $bred")
        end
    end

    # check if more jobs should be submitted
    for island in islands
        jobs_in_progress = round(Int,size(readdir("$(Rivanna.ROOT)/$island/prog"),1))
        println("In second loop $island jobs: $jobs_in_progress")
        if jobs_in_progress < MIN_WORK
            # below threshold!
            new_jobs = MIN_WORK - jobs_in_progress
            println("Making $new_jobs for $island")
            for j=1:new_jobs
                println("Now calling submitJob")
                Rivanna.submitJob(island)
                println("Submitted job")
            end
        end
    end

    # check if done population is too big
    for island in islands
        println("In third loop of do_maintainence")
        done_size = Rivanna.fetchDoneSize(island)
        println("In do_maintainence, $island has $(done_size) done jobs")
        if done_size > MAX_DONE
            println("$(done_size) > $(MAX_DONE)")
            gap = MAX_DONE - done_size
            deleted = Rivanna.cull(island)
            println("$island: culled $deleted")
        end
    end
end

println("GP SYSTEM FOR CIRIS\n")

if length(ARGS) == 1 && (ARGS[1] == "-d" || ARGS[1] == "--debug")
    # debug/default
    println("Using default arguments")
    root = "/scratch/cns7ae/archipelago"

    islands = ["G0"]
    Rivanna.initialize(root, islands)

    println("ROOT: $(Rivanna.ROOT)")
    println("ISLANDS: $(Rivanna.ISLANDS)")
elseif length(ARGS) >= 1
    Rivanna.initialize(ARGS[1], ARGS[2:end])
    println("ROOT: $(Rivanna.ROOT)")
    println("ISLANDS: $(Rivanna.ISLANDS)")
else
    # bad
    println("Bad call.")
    println("<path_to_root> island0 island1 island2 ...")
    exit(-1)
end

# the main loop
while true
    println("Starting loop again...")
    for i=1:5
        do_maintainence(Rivanna.ISLANDS)
        # sleep
        sleep(1)
        island_num = i - 1
    end
    println("Finishing loop")
end
