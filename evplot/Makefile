##-------------------------------------------------------------------##
##     Makefile for compiling with PGPLOT     			     ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2002	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FC =  g95
FFLAGS =  -Wall
LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng \
           -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation \
           -L/sw/lib/pgplot -lpgplot -lg2c

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $< -o $@

SOURCES= evsupersph.f90

OBJECTS = prompting.o transforms.o danpgtile.o system.o $(SOURCES:.f90=.o)

../evsupersph: $(OBJECTS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

danpgtile.o :
	$(FC) -c $(FFLAGS) ../plot/danpgtile.f -o $@ 

prompting.o :
	$(FC) -c $(FFLAGS) ../plot/prompting.f90 -o $@ 

transforms.o :
	$(FC) -c $(FFLAGS) ../plot/transform.f90 -o $@ 

system.o :
	$(FC) -c $(FFLAGS) ../plot/system_f2003.f90 -o $@ 
	
clean: 
	rm *.o *.mod
