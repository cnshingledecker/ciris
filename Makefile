# Compiler
#FC = gfortran
FC = ifort

# Flags
#FCFLAGS = -g  -march=native -fbounds-check -Wall -pg
#FCFLAGS =  -O3 -static-intel
#FCFLAGS = -O4 -march=native
FCFLAGS = -g -pg

OBJECTS = readinput.o subroutines.o functiondefs.o typedefs.o parameters.o main.o specdata.o gp.o branchmod.o bsimple.o

PROGRAM = ciris


$(PROGRAM): readinput.o subroutines.o functiondefs.o typedefs.o parameters.o main.o branchmod.o bsimple.o
	$(FC) -o $(PROGRAM) *.o $(FCFLAGS)

main.o: main.f90 subroutines.o parameters.o typedefs.o functiondefs.o gp.o branchmod.o bsimple.o readinput.o
	$(FC) -c main.f90 $(FCFLAGS)

typedefs.o: typedefs.f90
	$(FC) -c typedefs.f90 $(FCFLAGS)

specdata.o: specdata.f90 typedefs.o
	$(FC) -c specdata.f90 $(FCFLAGS)

readinput.o: readinput.f90 subroutines.o parameters.o functiondefs.o
	$(FC) -c readinput.f90 $(FCFLAGS)

subroutines.o: subroutines.f90 mc_toolbox.o typedefs.o functiondefs.o parameters.o specdata.o branchmod.o bsimple.o
	$(FC) -c subroutines.f90 $(FCFLAGS)

parameters.o: parameters.f90 typedefs.o
	$(FC) -c parameters.f90 $(FCFLAGS)

bsimple.o: bsimple.f90 parameters.o typedefs.o
	$(FC) -c bsimple.f90 $(FCFLAGS)

branchmod.o: branchmod.f90 parameters.o typedefs.o functiondefs.o
	$(FC) -c branchmod.f90 $(FCFLAGS)

mc_toolbox.o: mc_toolbox.f90
	$(FC) -c mc_toolbox.f90 $(FCFLAGS)

functiondefs.o: functiondefs.f90 parameters.o typedefs.o
	$(FC) -c functiondefs.f90 $(FCFLAGS)

gp.o: gp.f90 parameters.o typedefs.o subroutines.f90
	$(FC) -c gp.f90 $(FCFLAGS)

static: readinput.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o  $(FCFLAGS)

clean:
	rm -f *.mod *.o a.out
	rm -f fitness_results counter_test_wait_list.txt counter_test_wrong_spaces.txt
