#ifort -o  blastoff  chaco_data.f90 qbert.f90 moonbase.f90 parameters.f90  -march=native -O3 #0 -warn all -C 
gfortran -o  blastoff  mc_toolbox.f90 chaco_data.f03 qbert.f90 moonbase.f90  parameters.f90 functiondefs.f03 typedefs.f03  -march=native -O3 #0 -g  -pg -Wall -fbounds-check
