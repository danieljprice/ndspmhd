##-------------------------------------------------------------------##
##     Makefile for multirun.f90 	     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

F90C =  gfortran
F90FLAGS =  -O3 -fdefault-real-8
# PGPLOT libraries
##LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot -lpgplot -lg2c

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

SOURCES= multirun_dtconvergence.f90
 
OBJECTS = dimen1D.o variablesND.o defaults.o eos.o \
          kernel_erbs_utils.o kernelND.o readwrite_infile.o \
	  $(SOURCES:.f90=.o)

multirun: $(OBJECTS)	
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

dimen1D.o :
	$(F90C) -c $(F90FLAGS) ../src/dimen1D.f90 -o $@

eos.o :
	$(F90C) -c $(F90FLAGS) ../src/eos.f90 -o $@

kernelND.o :
	$(F90C) -c $(F90FLAGS) ../src/kernelND.f90 -o $@

kernel_erbs_utils.o :
	$(F90C) -c $(F90FLAGS) ../src/kernel_erbs_utils.f90 -o $@

variablesND.o :
	$(F90C) -c $(F90FLAGS) ../src/variablesND.f90 -o $@

defaults.o :
	$(F90C) -c $(F90FLAGS) ../src/defaults.f90 -o $@

readwrite_infile.o :
	$(F90C) -c $(F90FLAGS) ../src/readwrite_infile.f90 -o $@

clean:
	rm *.o *.mod
