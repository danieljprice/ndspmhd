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
LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11 -lF77

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

SOURCES= evsupersph.f90

OBJECTS = prompting.o danpgtile.o int_from_string.o $(SOURCES:.f90=.o)

../evsupersph: $(OBJECTS)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

danpgtile.o :
	$(F90C) -c $(F90FLAGS) ../plot/danpgtile.f -o $@ 

prompting.o :
	$(F90C) -c $(F90FLAGS) ../plot/prompting.f90 -o $@ 

int_from_string.o :
	$(F90C) -c $(F90FLAGS) ../plot/int_from_string.f90 -o $@ 
	
clean: 
	rm *.o
