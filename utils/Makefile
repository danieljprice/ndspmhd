##-------------------------------------------------------------------##
##     Makefile for compiling kernel plotting utilities	             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2004	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f95
F90C = f95

F90FLAGS = -C -O4
LDFLAGS = -lpgplot -L/usr/X11R6/lib -lX11 -lg2c -lpng ##-lF77
##LDFLAGS = -lpgplot -lX11 -lF77

# define the implicit rule to make a .o file from a .f90 file
%.o : %.f90
	$(F90C) -c $(F90FLAGS) $< -o $@

# define the implicit rule for compiling directly from a single .f90 file
##%.f90 :
##	$(F90C) -c $(F90FLAGS) $< -o $@

OBJECTS = dimen1D.o variablesND.o kernelND.o  \
	  kernelplot.o kernelstability1D.o danpgtile.o

kernelplot : $(OBJECTS)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

dimen1D.o :
	$(F90C) -c $(F90FLAGS) ../src/dimen1D.f90 -o $@

variablesND.o :
	$(F90C) -c $(F90FLAGS) ../src/variablesND.f90 -o $@

kernelND.o :
	$(F90C) -c $(F90FLAGS) ../src/kernelND.f90 -o $@

danpgtile.o :
	$(F90C) -c $(F90FLAGS) ../plot/danpgtile.f -o $@

clean :
	rm *.o *.mod
