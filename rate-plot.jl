println("Welcome to rate-plot.jl, where rates for CIRISS are plotted")
using DataFrames

path = "/home/cns/new_losalamos/"

reactions_outpule_filename = "ozone_reactions.csv"
reactionOutput = path*reactions_outpule_filename
println("The reactions output file is $reactionOutput")

speciesInput_filename     = "species.dat"
speciesInput = path*speciesInput_filename
println("The species input file is $speciesInput")

reactionsInput_filename   = "reactions.dat"
reactionsInput = path*reactionsInput_filename
println("The reactions input file is $reactionsInput")

#=
This function finds the index of something in a vector array in which it is
unique
=#
function findindex(p,itr)
  if p == "0"
    return 0
  end
  for n = 1:size(itr,1)
    if p == itr[n]
      return n
    end
  end
  return 666
end

#=
This functions Counts how many of something is in a vector array
=#
function howmany(thing,list)
  count = 0
  for n = 1:size(list,1)
    if thing == list[n]
      count = count + 1
    end
  end
  return count
end

function in_check(in_r,in_p,out_r,out_p)
  if size(union(in_r,out_r),1) == size(unique(in_r),1)
    test1 = true
  else
    test1 = false
  end

  if size(union(in_p,out_p),1) == size(unique(in_p),1)
    test2 = true
  else
    test2 = false
  end

  return test1*test2
end


#=
This function determines the class of a reactions
=#
function reactionclass(reactants,products,species)
  testReact = howmany(species,reactants)
#  println("testReact=$testReact")
  testProd = howmany(species,products)
#  println("testProd=$testProd")
  if testReact > testProd
    return "Destruction"
  elseif testProd > testReact
    return "Production"
  else
    return "Neutral"
  end
end

#=
Generate the dataframes containing the input and output dataframes
=#
println("Generating initial dataframes")
ozoneReactions = readtable(
                            reactionOutput,
                            separator = ',',
                            header = true,
#                            names = [:R1,:R2,:P1,:P2,:P3,:Fluence],
                            eltypes = [Int64,Int64,Int64,Int64,Int64,Float64]
                            )
println("Generated dataframe of reactions output")

species = readtable(
                    speciesInput,
                    names = [:Name, :E_D],
                    eltypes = [UTF8String,Float64],
                    header = false
)
println("Generated species input dataframe")

reactions = readtable(
                      reactionsInput,
                      header = false,
                      names = [:R1,:R2,:P1,:P2,:P3],
                      eltypes = [UTF8String,UTF8String,UTF8String,UTF8String,UTF8String]
)
println("Generated reactions input dataframe")

#=
Add reactions Type column to reactions dataframe and determine which one
each Type is
=#
println("Adding new column to reactions dataframe")
reactions[:Type] = "Neutral"
for i in 1:size(reactions[:Type],1)
  reactants = [reactions[:R1][i],reactions[:R2][i]]
  products  = [reactions[:P1][i],reactions[:P2][i],reactions[:P3][i]]
  reactions[:Type][i] = reactionclass(reactants,products,"O3")
end


println("Adding Count column to reactions dataframe")
reactions[:Count] = 0
fluenceLimit = 1.0E5
fluenceArr = zeros(0)
rateAnalytics = Array(Int64,(0,size(reactions[:R1],1)))
#=
Go through the reactions output, add any new reactions, and Count up rates
=#
for i in 1:size(ozoneReactions[:R1],1)
  outReacts = [ozoneReactions[:R1][i];ozoneReactions[:R2][i]]
  outProds  = [ozoneReactions[:P1][i];ozoneReactions[:P2][i];ozoneReactions[:P3][i]]
  in_network = false
  # Go through reactions for each line in the output reaction file
  for j in 1:size(reactions[:R1],1)
    in_reacts = [findindex(reactions[:R1][j],species[:Name]);
                 findindex(reactions[:R2][j],species[:Name])]
    inProds  = [findindex(reactions[:P1][j],species[:Name]);
                 findindex(reactions[:P2][j],species[:Name]);
                 findindex(reactions[:P3][j],species[:Name])]
    in_network = in_check(in_reacts,inProds,outReacts,outProds)
    if in_network == true
      reactions[:Count][j] += 1
      break
    end
  end
  # If the reaction isn't in the network, add it
  if in_network == false
    println("$outReacts, $outProds not in network")
    r1 = species[:Name][outReacts[1]]
    r2 = species[:Name][outReacts[2]]
    if outProds[1] == 0
      p1 = "0"
    else
      p1 = species[:Name][outProds[1]]
    end
    if outProds[2] == 0
      p2 = "0"
    else
      p2 = species[:Name][outProds[2]]
    end
    if outProds[3] == 0
      p3 = "0"
    else
      p3 = species[:Name][outProds[3]]
    end
    # Update the reactions dataframe
    push!(reactions,[r1;r2;p1;p2;p3;reactionclass(outReacts,outProds,7);0])
    # Expand the rateAnalytics array
    rateAnalytics = hcat(rateAnalytics,zeros(Int64,size(rateAnalytics,1)))
  end
  # Test to see if the fluence condition has been met, if so, save outputs and
  # re-initialize Counters
  if ozoneReactions[:Fluence][i] > fluenceLimit
    rateAnalytics = vcat(rateAnalytics,transpose(reactions[:Count]))
    reactions[:Count] = 0
    fluenceLimit += 0.1*fluenceLimit
    fluenceArr = vcat(fluenceArr,fluenceLimit)
    println("New fluence limit = $fluenceLimit")
  end
end

# Get the indices of all production and destruction reactions
prodReacts = Array(Int64,0)
destReacts = Array(Int64,0)
neutReacts = Array(Int64,0)
for i in 1:size(reactions[:Type],1)
  if reactions[:Type][i] == "Production"
    push!(prodReacts,i)
  elseif reactions[:Type][i] == "Destruction"
    push!(destReacts,i)
  elseif reactions[:Type][i] == "Neutral"
    push!(neutReacts,i)
  end
end
println("prodReacts=$prodReacts")
println("destReacts=$destReacts")
println("neutReacts=$neutReacts")

println("Now initializing the percentRates array")
totProdArr = zeros(0)
totDestArr = zeros(0)
totArr     = zeros(0)
percentRates = zeros(Float64,(size(rateAnalytics,1),size(rateAnalytics,2)))
println("percentRates array now initialized")
for i in 1:size(rateAnalytics,1)
  totalProd = 0
  totalDest = 0
  totalNeut = 0
  # Go through row and add up reaction types
  for j in 1:size(rateAnalytics,2)
    if in(j,prodReacts)
      totalProd += rateAnalytics[i,j]
    elseif in(j,destReacts)
      totalDest += rateAnalytics[i,j]
    elseif in(j,neutReacts)
      totalNeut += rateAnalytics[i,j]
    end
  end
  # Go back through row and get percentage rates
  for j in 1:size(rateAnalytics,2)
    if in(j,prodReacts)
      if totalProd != 0
        percentRates[i,j] = Float64(rateAnalytics[i,j])/Float64(totalProd)
      end
    elseif in(j,destReacts)
      if totalDest != 0
        percentRates[i,j] = Float64(rateAnalytics[i,j])/Float64(totalDest)
      end
    elseif in(j,neutReacts)
      if totalNeut != 0
        percentRates[i,j] = Float64(rateAnalytics[i,j])/Float64(totalNeut)
      end
    end
  end
  totProdArr = vcat(totProdArr,totalProd)
  totDestArr = vcat(totDestArr,totalDest)
  totArr     = vcat(totArr,totalProd+totalDest)
end

proddf = DataFrame()
destdf = DataFrame()
neutdf = DataFrame()
println("Plot dataframes initialized")
for i in 1:size(percentRates,2)
  r1 = reactions[:R1][i]
  r2 = reactions[:R2][i]
  prods = [reactions[:P1][i]; reactions[:P2][i];reactions[:P3][i]]
  rname = "$r1 + $r2 -> "
  pcount = 0
  for q in 1:3
    if pcount == 0
      if prods[q] != "0"
        rname = rname*"$(prods[q]) "
        pcount += 1
      end
    else
      if prods[q] != "0"
        rname = rname*"+ $(prods[q])"
      end
    end
  end
  if in(i,prodReacts)
    proddf = vcat(proddf,DataFrame(
                                   Fluence=fluenceArr,
                                   PercentTotal=percentRates[:,i],
                                   label="$rname"
                                   ))
  elseif in(i,destReacts)
    destdf = vcat(destdf,DataFrame(
                                   Fluence=fluenceArr,
                                   PercentTotal=percentRates[:,i],
                                   label="$rname"
                                   ))
  elseif in(i,neutReacts)
    neutdf = vcat(neutdf,DataFrame(
                                   Fluence=fluenceArr,
                                   PercentTotal=percentRates[:,i],
                                   label="$rname"
                                   ))
  end
end

println("Now generating plots")
p1 = plot(
          proddf[1e16.>proddf[:Fluence].>1e12,:],
          x="Fluence",
          y="PercentTotal",
          color="label",
          Geom.smooth,
          Guide.title("Production Reactions"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("Fraction of total"),
#          Coord.Cartesian(xmin=fluenceArr[1],xmax=fluenceArr[size(fluenceArr,1)])
          Scale.x_log10(minvalue=1E12,maxvalue=1e16)
         )


p2 = plot(
          destdf[1e16.>destdf[:Fluence].>1e12,:],
          x="Fluence",
          y="PercentTotal",
          color="label",
          Geom.smooth,
          Guide.title("Destruction Reactions"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("Fraction of total"),
#          Coord.Cartesian(xmin=fluenceArr[1],xmax=fluenceArr[size(fluenceArr,1)])
          Scale.x_log10(minvalue=1E12,maxvalue=1e16)
         )



totProdPercent = zeros(0)
totDestPercent = zeros(0)

for i in 1:size(totProdArr,1)
  if totArr[i] != 0
    totProdPercent = vcat(totProdPercent,totProdArr[i]/totArr[i])
    totDestPercent = vcat(totDestPercent,totDestArr[i]/totArr[i])
  else
    totProdPercent = vcat(totProdPercent,0.0)
    totDestPercent = vcat(totDestPercent,0.0)
  end
end
totdf = vcat(
             DataFrame(
                       Fluence=fluenceArr,
                       PercentTotal=totProdPercent,
                       label="Total production"
             ),
             DataFrame(
                       Fluence=fluenceArr,
                       PercentTotal=totDestPercent,
                       label="Total destruction"
             )
)


p3 = plot(
          totdf[1e16.>totdf[:Fluence].>1e12,:],
          x="Fluence",
          y="PercentTotal",
          color="label",
          Geom.smooth,
          Guide.title("Relative Production and Destruction"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("Fraction of total"),
          #Coord.Cartesian(xmin=fluenceArr[1],xmax=fluenceArr[size(fluenceArr,1)])
          Scale.x_log10(minvalue=1E12,maxvalue=1e16)
         )

println("Generating final output pdf")
draw(PDF("/home/cns/new_losalamos/rate_plots.pdf",6inch,9inch),vstack(p1,p2,p3))

println("Ending rate-plot.jl!")
