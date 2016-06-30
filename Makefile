# Compiler
  FC = gfortran
  #FC = ifort

# Flags
  #FCFLAGS = -g  -march=native -fbounds-check -Wall
  #FCFLAGS =  -O3 -static-intel 
  FCFLAGS = -O3 -march=native

OBJECTS = qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o specdata.o gp.o

PROGRAM = losalamos


$(PROGRAM): qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o $(FCFLAGS)

main.o: main.f90 subroutines.o parameters.o typedefs.o functiondefs.o gp.o
	$(FC) -c main.f90

typedefs.o: typedefs.f90
	$(FC) -c typedefs.f90

specdata.o: specdata.f90 typedefs.o
	$(FC) -c specdata.f90

qbert.o: qbert.f90 subroutines.o parameters.o
	$(FC) -c qbert.f90

subroutines.o: subroutines.f90 mc_toolbox.o typedefs.o functiondefs.o parameters.o specdata.o
	$(FC) -c subroutines.f90

parameters.o: parameters.f90
	$(FC) -c parameters.f90

mc_toolbox.o: mc_toolbox.f90
	$(FC) -c mc_toolbox.f90

functiondefs.o: functiondefs.f90 parameters.o typedefs.o
	$(FC) -c functiondefs.f90

gp.o: gp.f90 parameters.o typedefs.o
	$(FC) -c gp.f90

static: qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o -O3 -static-intel 

clean:
	rm -f *.mod *.o a.out
	rm -f fitness_results counter_test_wait_list.txt counter_test_wrong_spaces.txt
