##-------------------------------------------------------------------##
##     Makefile for compiling with PGPLOT     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

F90C =  g95
F90FLAGS =  -Wall
LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng \
           -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation \
           -L/sw/lib/pgplot -lpgplot -lg2c

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@

SOURCES= evsupersph.f90

OBJECTS = prompting.o transforms.o danpgtile.o system.o $(SOURCES:.f90=.o)

../evsupersph: $(OBJECTS)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

danpgtile.o :
	$(F90C) -c $(F90FLAGS) ../plot/danpgtile.f -o $@ 

prompting.o :
	$(F90C) -c $(F90FLAGS) ../plot/prompting.f90 -o $@ 

transforms.o :
	$(F90C) -c $(F90FLAGS) ../plot/transform.f90 -o $@ 

system.o :
	$(F90C) -c $(F90FLAGS) ../plot/system_f2003.f90 -o $@ 
	
clean: 
	rm *.o *.mod
