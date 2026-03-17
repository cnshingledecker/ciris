# Compiler
FC = gfortran
#FC = ifort

# Flags
#FCFLAGS = -g -march=native -fbounds-check -Wall
FCFLAGS = -march=native -O3 -ffast-math

PROGRAM = losalamos

OBJECTS = typedefs.o parameters.o mc_toolbox.o functiondefs.o specdata.o \
          gp.o subroutines.o qbert.o bresenham.o main.o

.PHONY: all static clean

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS)

static: $(OBJECTS)
	$(FC) $(FCFLAGS) -static -o $(PROGRAM) $(OBJECTS)

# Pattern rule — covers all .f03 -> .o compilations
%.o: %.f03
	$(FC) $(FCFLAGS) -c $<

# Explicit dependency ordering (for Fortran module files)
main.o:         main.f03        subroutines.o parameters.o typedefs.o functiondefs.o gp.o
typedefs.o:     typedefs.f03
specdata.o:     specdata.f03    typedefs.o
qbert.o:        qbert.f03       subroutines.o parameters.o
subroutines.o:  subroutines.f03 mc_toolbox.o typedefs.o functiondefs.o parameters.o specdata.o
parameters.o:   parameters.f03
mc_toolbox.o:   mc_toolbox.f03
functiondefs.o: functiondefs.f03 parameters.o typedefs.o
gp.o:           gp.f03          parameters.o
bresenham.o:    bresenham.f03

clean:
	rm -f *.mod *.o $(PROGRAM) a.out
	rm -f fitness_results counter_test_wait_list.txt counter_test_wrong_spaces.txt
