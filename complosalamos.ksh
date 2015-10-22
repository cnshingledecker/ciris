#ifort -o  losalamos  chaco_data.f90 qbert.f90 moonbase.f90 parameters.f90  -march=native -O3 #0 -warn all -C 
gfortran -o  losalamos mc_toolbox.f03 chaco_data.f03 qbert.f03 moonbase.f03  parameters.f03 functiondefs.f03 typedefs.f03  -march=native -O0 -g  -pg -Wall -fbounds-check -static-libgcc
