print("Welcome to abplot, where abundances produced by CIRISS are plotted...")

filename = '/home/cns/losalamos/abundance.csv'

ciriss_output = read.csv(filename, 
                         header = FALSE, 
                         sep = ",", 
                         col.names = c("Time", "Fluence", "O", "O3", "Protons", "O3 prod/dest"), 
                         colClasses = c(rep("numeric",6))
                         )


plot.default(ciriss_output$Fluence,ciriss_output$O3,
     xlab=expression("Fluence (" ~ H^{phantom()+phantom()} ~ cm^{phantom()-2} ~ ")" ),
     ylab=expression("[O"[3] ~ "]"),
     main=expression("[O"[3] ~ "] abundance vs Fluence"),
     log = "x"
     )

par(mar = c(5.1, 8.1, 5.1, 8.1))