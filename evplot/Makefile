##-------------------------------------------------------------------##
##     Makefile for compiling with PGPLOT     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
FC = f90

FFLAGS = 
LDFLAGS = -lpgplot -L/local/lib/pgplot -lX11 -lF77

SOURCES= evsupersph.f 

OBJECTS = $(SOURCES:.f=.o)

../evsupersph: $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

