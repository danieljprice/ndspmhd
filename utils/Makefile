##-------------------------------------------------------------------##
##     Standard Makefile for f90        			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
F90C = f90

F90FLAGS = 
LDFLAGS = -lpgplot -lX11 -lF77
# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

SOURCES= sametime.f90

OBJECTS = $(SOURCES:.f90=.o)

OBJECTS2 = ../src/dimen1D.o ../src/variablesND.o ../src/kernelND.o kernelplot.o

sametime : $(OBJECTS)
	$(F90C) $(F90FLAGS) -o $@ $(OBJECTS)

kernelplot : $(OBJECTS2)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS2)
