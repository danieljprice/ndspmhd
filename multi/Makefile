##-------------------------------------------------------------------##
##     Makefile for multirun.f90 	     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
F90C = f90

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

SOURCES= multirun.f90
 
OBJECTS = dimen1D.o variablesND.o defaults.o \
          read_infile_mhd.o write_infile_mhd.o int_from_string.o \
	  $(SOURCES:.f90=.o)

multirun: $(OBJECTS)	
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

dimen1D.o :
	$(F90C) -c $(F90FLAGS) ../src/dimen1D.f90 -o $@

variablesND.o :
	$(F90C) -c $(F90FLAGS) ../src/variablesND.f90 -o $@

defaults.o :
	$(F90C) -c $(F90FLAGS) ../src/defaults.f90 -o $@

read_infile_mhd.o :
	$(F90C) -c $(F90FLAGS) ../src/read_infile_mhd.f90 -o $@

write_infile_mhd.o :
	$(F90C) -c $(F90FLAGS) ../src/write_infile_mhd.f90 -o $@

int_from_string.o :
	$(F90C) -c $(F90FLAGS) ../plot/int_from_string.f90 -o $@

clean:
	rm *.o *.mod
