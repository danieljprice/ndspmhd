##-------------------------------------------------------------------##
##     Makefile for multirun.f90 	     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##
##    First time this is run, type Make initial
##    to compile the relevant files from the SPMHD code
##    After that, just make as usual
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
F90C = lf95 	##f90

# use -C for array out-of-bounds checking
# -u forces implicit none
FFLAGS = -04 -cg92 -lm
#F90FLAGS = -O4 -cg92 -lm 
# use the following instead for double precision
F90FLAGS = -O	##4 -cg92 -lm -C -xtypemap=real:64,double:64,integer:32
# PGPLOT libraries
#LDFLAGS = -lpgplot -L/local/lib/pgplot -lX11 -lF77


# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@


SOURCES= multirun.f90 int_from_string.f90

SOURCESCODE = ../src/dimen1D.f90\
	 ../src/variablesND.f90 \
	 ../src/read_infile_mhd.f90 \
	 ../src/write_infile_mhd.f90 	 
 
OBJECTS = $(SOURCES:.f90=.o) dimen1D.o variablesND.o read_infile_mhd.o write_infile_mhd.o

multirun: $(OBJECTS)	
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

initial: $(SOURCESCODE)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -c $(SOURCESCODE)

clean:
	rm *.o *.mod
