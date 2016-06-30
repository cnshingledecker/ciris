using DataFrames
using Gadfly

# Input parameters
path = "/home/cns/new_losalamos/"
abundanceFile = "abundance.csv"
plotname = "abundance_plot.pdf"
thick = 3.0e-6
edge  = 3.0e-6
denom = thick*edge*edge
iceSize = (102^3)/3

println("Reading input file")
modelOutput = readtable(
                  path*abundanceFile,
                  names = [:Time, :Fluence, :O, :O3, :ProtonNum, :ProdDest],
                  eltypes = [Float64,Float64,Float64,Float64,Int64,Float64],
                  header = false
                  )

println("Initializing dataframes")
O3ab_df = DataFrame()
Oab_df = DataFrame()
O3frac_df = DataFrame()
Ofrac_df = DataFrame()

println("Populating dataframes")
Oab_df = vcat(Oab_df,DataFrame(
                             Fluence=modelOutput[:Fluence],
                             Abundance=1.0e20*(modelOutput[:O]/denom),
                             Species="O")
            )


Ofrac_df = vcat(Ofrac_df,DataFrame(Fluence=modelOutput[:Fluence],
                         Abundance=modelOutput[:O]/iceSize,
                         Species="O")
               )


O3ab_df = vcat(O3ab_df,DataFrame(
                             Fluence=modelOutput[:Fluence],
                             Abundance=1.0e20*(modelOutput[:03]/denom),
                             Species="O3")
               )


O3frac_df = vcat(O3frac_df,DataFrame(
                             Fluence=modelOutput[:Fluence],
                             Abundance=modelOutput[:03]/iceSize,
                             Species="O3")
                 )
println("Generating plots")
p1 = plot(
          O3ab_df[1e16.>O3ab_df[:Fluence].>1e12,:],
          x="Fluence",
          y="Abundance",
          color="Species",
          Geom.smooth,
          Guide.title("[O3] vs. Fluence"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("[O3]*1e20"),
          Scale.x_log10(minvalue = 1e12, maxvalue = 1e16)
)

p2 = plot(
          O3frac_df[1e16.>O3frac_df[:Fluence].>1e12,:],
          x="Fluence",
          y="Abundance",
          color="Species",
          Geom.smooth,
          Guide.title("O3 Fraction of total ice vs. Fluence"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("Fraction"),
          Scale.x_log10(minvalue = 1e12, maxvalue = 1e16)
)

p3 = plot(
          Oab_df[1e16.>Oab_df[:Fluence].>1e12,:],
          x="Fluence",
          y="Abundance",
          color="Species",
          Geom.smooth,
          Guide.title("[O] vs. Fluence"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("[O]*1e20"),
          Scale.x_log10(minvalue = 1e12, maxvalue = 1e16)
)

p4 = plot(
          Ofrac_df[1e16.>Ofrac_df[:Fluence].>1e12,:],
          x="Fluence",
          y="Abundance",
          color="Species",
          Geom.smooth,
          Guide.title("O Fraction of total ice vs. Fluence"),
          Guide.xlabel("Fluence"),
          Guide.ylabel("Fraction"),
          Scale.x_log10(minvalue = 1e12, maxvalue = 1e16)
)

println("Now exporting plot")
draw(PDF(path*plotname,6inch,9inch),vstack(p1,p2,p3,p4))
println("Ending script!!")
