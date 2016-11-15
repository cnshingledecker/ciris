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
MIN_WORK      = 10  # max number of jobs/island
MIN_TODO      = 200 # min todo size/island
MAX_DONE      = 30 # max done size/island
EXILE         = 25 # not implemented yet...
MUTATE_CHANCE = 3  # 1 out of...
EXILE_CHANCE  = 4  # 1 out of...

# Some sanity assertions
@assert MIN_WORK < MIN_TODO "Max possible jobs/island exceeds minimum population threshold!"

# control logic for maintaining jobs and islands
function do_maintainence(islands)
    # check if todo population should be bolstered
    temp_best_fitness = 9.9E12
    for island in islands
        todo_size = Rivanna.fetchTodoSize(island)
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
            cd("$(Rivanna.ROOT)/$island/done")
            for j=1:deficit
                if rand(1:MUTATE_CHANCE) == 1
                    #mutate
                    # pick a random file and fetch its contents
                    file_contents = "$(Rivanna.ROOT)/$island/done/$(rand(done_files))"
                    parsed_file = ParameterIO.readSolution(file_contents)
                    # drop fitness and seed
                    delete!(parsed_file, "FITNESS")
                    delete!(parsed_file, "SEED")
                    # vary number of mutations
                    for i = rand(1:4)
                        ParameterIO.mutate!(parsed_file)
                    end
                    name = ParameterIO.writeSolution(parsed_file)
                    run(`cp $(Rivanna.ROOT)/$island/done/$name $(Rivanna.ROOT)/$island/todo`)
                    run(`rm $(Rivanna.ROOT)/$island/done/$name`)
                    mutated = mutated + 1
                else
                    #breed
                    # pick two files at random
                    f1, f2 = rand(done_files, 2)
#                    f1_contents = Rivanna.fetchCat("$(Rivanna.ROOT)/$island/done/$f1")
                    f1_contents = "$(Rivanna.ROOT)/$island/done/$f1"
#                    f2_contents = Rivanna.fetchCat("$(Rivanna.ROOT)/$island/done/$f2")
                    f2_contents = "$(Rivanna.ROOT)/$island/done/$f2"
                    f1_parsed = ParameterIO.readSolution(f1_contents)
                    f2_parsed = ParameterIO.readSolution(f2_contents)
                    delete!(f1_parsed, "FITNESS")
                    delete!(f2_parsed, "FITNESS")
                    delete!(f1_parsed, "SEED")
                    delete!(f2_parsed, "SEED")
                    # breed and stash new candidate
                    child = ParameterIO.breed(f1_parsed, f2_parsed)
                    name = ParameterIO.writeSolution(child)
                    run(`cp $(Rivanna.ROOT)/$island/done/$name $(Rivanna.ROOT)/$island/todo`)
                    run(`rm $(Rivanna.ROOT)/$island/done/$name`)
                    bred = bred + 1
                end
            end
            println("$island: mutated $mutated, bred $bred")
        end
    end

    # check if more jobs should be submitted
    # Check to see if one needs to clean out the prog folder
    for island in islands
        jobs_in_progress = round(Int,size(readdir("$(Rivanna.ROOT)/$island/prog"),1))
        if jobs_in_progress < MIN_WORK
            # below threshold!
            new_jobs = MIN_WORK - jobs_in_progress
            for j=1:new_jobs
                Rivanna.submitJob(island)
            end
        end
    end

    # check if done population is too big
    for island in islands
        done_size = Rivanna.fetchDoneSize(island)
        if done_size > MAX_DONE
            gap = MAX_DONE - done_size
            deleted, temp_best_fitness = Rivanna.cull(island)
            println("$island: culled $deleted")
        end
    end
    return temp_best_fitness
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
BEST_FITNESS  = 9.9E12
temp_best = BEST_FITNESS
count_count = 0
while true
    println("************************************")
    println("BEST_FITNESS=$(BEST_FITNESS)")
    println("************************************")
        temp_best = do_maintainence(Rivanna.ISLANDS)
        for island in Rivanna.ISLANDS
          njobs = parse(Int,readstring(pipeline(`squeue -u cns7ae`,`wc -l`)))
          nprog = parse(Int,readstring(pipeline(`ls $(Rivanna.ROOT)/$island/prog`,`wc -l`)))
          # If the number of jobs is less than the maxjob, rm all jobs
          if njobs < nprog  
             println("$njobs in queue < $nprog in prog")
             run(`echo Deleting jobs in $(Rivanna.ROOT)/G0/prog`) 
             run(`ls $(Rivanna.ROOT)/G0/prog`) 
             run(`rm -rf $(Rivanna.ROOT)/G0/prog`) 
             run(`mkdir $(Rivanna.ROOT)/G0/prog`) 
          njobs = parse(Int,readstring(pipeline(`squeue -u cns7ae`,`wc -l`)))
          println("After Deletions!: There are $njobs running currently on $(length(Rivanna.ISLANDS))")
          nprog = parse(Int,readstring(pipeline(`ls $(Rivanna.ROOT)/$island/prog`,`wc -l`)))
          println("After Deletions:! There are $nprog current jobs in $island/prog")

          end
        end 
        if temp_best < BEST_FITNESS
            BEST_FITNESS = temp_best
        end
        # sleep
        sleep(30)
    println("Finishing loop")
end
