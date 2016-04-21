# Compiler
  FC = gfortran
  #FC = ifort

# Flags
  #FCFLAGS = -g  -march=native -fbounds-check -Wall
  FCFLAGS = -march=native -O3 -ffast-math

OBJECTS = qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o specdata.o gp.o

PROGRAM = losalamos


$(PROGRAM): qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o $(FCFLAGS)

main.o: main.f03 subroutines.o parameters.o typedefs.o functiondefs.o gp.o
	$(FC) -c main.f03

typedefs.o: typedefs.f03
	$(FC) -c typedefs.f03

specdata.o: specdata.f03 typedefs.o
	$(FC) -c specdata.f03

qbert.o: qbert.f03 subroutines.o parameters.o
	$(FC) -c qbert.f03

subroutines.o: subroutines.f03 mc_toolbox.o typedefs.o functiondefs.o parameters.o specdata.o
	$(FC) -c subroutines.f03

parameters.o: parameters.f03
	$(FC) -c parameters.f03

mc_toolbox.o: mc_toolbox.f03
	$(FC) -c mc_toolbox.f03

functiondefs.o: functiondefs.f03 parameters.o typedefs.o
	$(FC) -c functiondefs.f03

gp.o: gp.f03 parameters.o
	$(FC) -c gp.f03

static: qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o -O3 -ffast-math -static

clean:
	rm -f *.mod *.o a.out
	rm -f fitness_results counter_test_wait_list.txt counter_test_wrong_spaces.txt
