##-------------------------------------------------------------------##
##     Standard Makefile for f90        			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
FC = f90

FFLAGS = 
LDFLAGS = 

SOURCES= 

OBJECTS = $(SOURCES:.f=.o)

.f.o : $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

