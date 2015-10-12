#ifort -o  blastoff  chaco_data.f90 qbert.f90 moonbase.f90 parameters.f90  -march=native -O0 -pg #0 -C -g  -traceback

gfortran -o  transport  chaco_data.f90 qbert.f90 transport_test.f90  parameters.f90   -march=native -O3 #0 -g  -Wall -fbounds-check
