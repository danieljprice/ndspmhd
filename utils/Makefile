##-------------------------------------------------------------------##
##     Makefile for compiling kernel plotting utilities	             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2004	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler flags
FORTRAN = ifc
F90C =  ifc
FC = ifc
F90FLAGS =  -O -WB ##-C
LDFLAGS =  -lpgplot -L/usr/X11R6/lib -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2/ -lg2c -lpng -Vaxlib

.SUFFIXES:
.SUFFIXES:.f90 .o
# Fortran flags same as F90
#FC = $(F90C)
#FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file
%.o : %.f90
	$(F90C) -c $(F90FLAGS) $< -o $@

# define the implicit rule for compiling directly from a single .f90 file
#%.f90 :
#	$(F90C) -c $(F90FLAGS) $< -o $@

OBJECTS = dimen1D.o variablesND.o kernelND.o  \
	  kernelplot.o kernelstability1D.o danpgtile.o

OBJSPHDERIV = kernelND.o variablesND.o sphderiv.o interpolate_kernel.o

kernelplot : $(OBJECTS)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

sametime : sametime.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ sametime.o

sphderiv : $(OBJSPHDERIV)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJSPHDERIV)

dimen1D.o :
	$(F90C) -c $(F90FLAGS) ../src/dimen1D.f90 -o $@

variablesND.o :
	$(F90C) -c $(F90FLAGS) ../src/variablesND.f90 -o $@

kernelND.o :
	$(F90C) -c $(F90FLAGS) ../src/kernelND.f90 -o $@

danpgtile.o :
	$(F90C) -c $(F90FLAGS) ../plot/danpgtile.f -o $@

interpolate_kernel.o :
	$(F90C) -c $(F90FLAGS) ../src/interpolate_kernel.f90 -o $@


clean :
	rm *.o *.mod
