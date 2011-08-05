##-------------------------------------------------------------------##
##     Makefile for compiling with PGPLOT     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
F90C = f90

FFLAGS = 
F90FLAGS = -O
LDFLAGS = -lpgplot -lX11 -lF77

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

SOURCES= evsupersph.f90

OBJECTS = $(SOURCES:.f90=.o) danpgtile.o

../evsupersph: $(OBJECTS)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

danpgtile.o :
	$(F90C) -c $(F90FLAGS) ../plot/danpgtile.f -o $@ 
	
clean: 
	rm *.o
