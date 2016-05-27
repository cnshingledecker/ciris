library(scatterplot3d)

print("Welcome to Trackplot, where tracks are plotted...")
filename <- "/home/cns/losalamos/trackplot.csv"
range <- c(1,340)
track_coords <- read.csv(filename,
                         header = FALSE, 
                         sep = ",", 
                         col.names = c("z","y","x"),
                         colClasses = c(rep("integer",3))
                         )

scatterplot3d(x=track_coords$x,
              y=track_coords$y,
              z=track_coords$z,
              xlim=range,
              ylim=range,
              zlim=range,
              pch=("."),
              xlab = "x",
              ylab = "y",
              zlab = "z"
              )