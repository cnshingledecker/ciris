module ParameterIO

using Base.Random.uuid4

constraints = Dict(
    "O" => 55:0.01:400,               #1
    "O3" => 100:0.01:800,             #2
    "ELASTIC_LOSS" => 0.0:0.0001:0.5, #3
    "TRL_NU" => 1E10:1E10:1E12,       #4
    "O2_DISPROB" => 0.0:0.001:1.0,    #5
    "O3_DISPROB" => 0.0:0.001:1.0,    #6
    "STEPFAC" => 1E-3:1E-4:1.0,       #7
    "ESTEPFAC" => 1E-3:1E-3:1.0,      #8
    "AVAL" => 1:1:60,                 #9
)
parameters = collect(keys(constraints))

function sanitycheck(arg)
    println("sanitycheck with arg ", arg)
end

##################################
### Breed a Candidate Solution ###
##################################
#NOTE!  This methodology creates a bias towards the center of
#       the solution space
#       Maybe something like a probability distribution?
function breed(solution1, solution2)
    child = mixTraits(solution1, solution2)
    return child
end

# takes a random combination of traits from each parent
function mixTraits(solution1, solution2)
    child = copy(solution1)
    for i = keys(child)
        if rand(0:1) == 1
            child[i] = solution2[i]
        end
    end
    return child
end

# takes two values from parents and their container c
# returns a new trait
function breedValue(val1::AbstractFloat, val2::AbstractFloat, c::FloatRange)
    # if val is a subtype of Array aka set of things
    min = val1
    max = val2
    if val1 > val2
        tmp = val1
        val1 = val2
        val2 = tmp
    end
    return min + rand() * (max-min)
end

function breedValue(val1::Integer, val2::Integer, c::UnitRange)
    min = val1
    max = val2
    if min > max
        tmp = min
        min = max
        max = tmp
    end
    return rand(min:max)
end

###################################
### Mutate a Candidate Solution ###
###################################
function mutate(solution)
    new_sol = copy(solution)
    mutate!(new_sol)
    return new_sol
end

function mutate!(solution)
    random_attribute = rand(parameters)
    solution[random_attribute] = generateRandom(constraints, random_attribute)
    return solution
end

""" Convenience function for end to end mutation
    returns name of file
"""
function mutateFile(filename)
    path = join(split(filename, '/')[1:end-1], '/')
    s = readSolution(filename)
    mutate!(s)
    name = string(path, '/', string(uuid4()))
    writecsv(name, s)
    return name
end

#######################################
### Generate new Candidate Solution ###
#######################################
function generateSolution()
    sol = Dict()
    for (i,j)=ParameterIO.constraints
        if i == "FAST_REACTS" || i == "FRAGILE"
            sol[i] = """'$(generateRandom(constraints, i))'"""
        else
            sol[i] = generateRandom(constraints, i)
        end
    end
    return sol
end

""" Generates a random value using a range resolved by a dict and
    a key. A FloatRange defaults to step 1, work around...
"""
function generateRandom(dict, key)
    if dict[key] == 0:1
        rand()
    else
        rand(dict[key])
    end
end


#############################
### Candidate Solution IO ###
#############################
function readSolution(filename)
    file = readcsv(filename)
    return Dict(zip(file[:,1], file[:,2]))
end

function writeSolution(solution)
    name = string(uuid4())
    writecsv(name, solution)
    return name
end


#####################################
### Simulation Parameter Files IO ###
#####################################
function readSpeciesFile(filename)
    file = readcsv(filename)
    return Dict(zip(file[:,1], file[:,2]))
end

function writeSpeciesFile(filename, dict)
    writecsv(filename, dict)
end

end
