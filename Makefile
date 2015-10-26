# Compiler
  #FC = gfortran
  FC = gfortran5
  #FC = ifort

# Flags
  #FCFLAGS = -g -pg -static-libgcc -march=native -fbounds-check -Wall
  FCFLAGS = -march=native -O3 -ffast-math -static-libgcc

OBJECTS = qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o specdata.o

PROGRAM = losalamos


$(PROGRAM): qbert.o subroutines.o functiondefs.o typedefs.o parameters.o main.o
	$(FC) -o $(PROGRAM) *.o $(FCFLAGS)

main.o: main.f03 subroutines.o parameters.o typedefs.o functiondefs.o
	$(FC) -c main.f03

typedefs.o: typedefs.f03
	$(FC) -c typedefs.f03

#typedefs.mod: typedefs.f03 typedefs.o
#	@true

specdata.mod: specdata.f03 specdata.o
	@true


specdata.o: specdata.f03 typedefs.o
	$(FC) -c -M specdata.f03 

qbert.o: qbert.f03 subroutines.o parameters.o
	$(FC) -c qbert.f03

subroutines.o: subroutines.f03 mc_toolbox.o typedefs.o functiondefs.o parameters.o
	$(FC) -c subroutines.f03



mc_toolbox.o: mc_toolbox.f03
<<<<<<< HEAD
	$(FC) -c $< 

functiondefs.o: functiondefs.f03 parameters.o typedefs.o
	$(FC) -c $< 


parameters.o: parameters.f03 
	$(FC) -c $< 
=======
	$(FC) -c mc_toolbox.f03 

functiondefs.o: functiondefs.f03 parameters.o typedefs.o
	$(FC) -c functiondefs.f03 


parameters.o: parameters.f03 
	$(FC) -c parameters.f03 
>>>>>>> created 1.1

#parameters.mod: parameters.f03  parameters.o 
#	@true

clean:
	rm -f *.mod *.o
