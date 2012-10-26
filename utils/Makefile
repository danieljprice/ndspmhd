##-------------------------------------------------------------------##
##     Makefile for compiling various utilities which use	     ##
##      parts of the main code                                       ##
##                                                                   ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, -2004	 	     ##
##     University of Exeter, Exeter, UK       2004-                  ##
##-------------------------------------------------------------------##


.KEEP_STATE:

## Compiler flags
F90C =  gfortran
F90FLAGS =  -O3 -Wall -fbounds-check
#LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot95 -lpgplot
LDFLAGS =  -L/usr/X11R6/lib -lX11 -L${PGPLOT_DIR} -lpgplot -lpng
DBLFLAGS = -fdefault-real-8
ifndef SPLASH_DIR
SPLASH_DIR=$(HOME)/splash/src
endif

.SUFFIXES:
.SUFFIXES:.f .f90 .o .mod
VPATH=../src $(SPLASH_DIR)
# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)
CC = $(F90C)

# define the implicit rule to make a .o file from a .f90 file
%.o : %.f90
	$(F90C) -c $(F90FLAGS) $< -o $@

# define the implicit rule for compiling directly from a single .f90 file
#%.f90 :
#	$(F90C) $(F90FLAGS) $(LDFLAGS) $< -o $@

OBJECTS=setpage.o kernel_cspline_utils.o kernel_erbs_utils.o kernelND.o legend.o

OBJSPHDERIV = kernelND.o sphderiv.o

OBJMODDUMP = variablesND.o allocateND.o prompting.o \
             readwrite_dumps.o modifydump.o

all: kernelplot kernelplot3D getkernellist

kernelplot: $(OBJECTS) kernelplot.o kernelstability1D.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS) kernelplot.o kernelstability1D.o

kernelplot3D: $(OBJECTS) kernelplot3D.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS) kernelplot3D.o

getkernellist: $(OBJECTS) getkernellist.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJECTS) getkernellist.o

kernelplot.o : kernelND.o

sametime : sametime.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ sametime.o

plotfreq : plotfreq.o danpgtile.o
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ plotfreq.o danpgtile.o

plotfreq.o : plotfreq.f90
	$(F90C) -c $(F90FLAGS) plotfreq.f90 -o $@

sphderiv : $(OBJSPHDERIV)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o $@ $(OBJSPHDERIV)
        
moddump1D: dimen1D.o $(OBJMODDUMP)
	$(F90C) $(F90FLAGS) -o $@ dimen1D.o $(OBJMODDUMP)

moddump2D: dimen2D.o $(OBJMODDUMP)
	$(F90C) $(F90FLAGS) -o $@ dimen2D.o $(OBJMODDUMP)

moddump3D: dimen3D.o $(OBJMODDUMP)
	$(F90C) $(F90FLAGS) -o $@ dimen3D.o $(OBJMODDUMP)

dimen1D.o :
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/dimen1D.f90 -o $@

dimen2D.o :
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/dimen2D.f90 -o $@

dimen3D.o :
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/dimen3D.f90 -o $@

variablesND.o :
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/variablesND.f90 -o $@

kernelND.o :
	$(F90C) -c $(F90FLAGS) ../src/kernelND.f90 -o $@

allocateND.o : variablesND.o
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/allocateND.f90 -o $@

readwrite_dumps.o : variablesND.o geometry.o
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../src/readwrite_dumps.f90 -o $@

modifydump.o : modifydump.f90
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) modifydump.f90 -o $@

geometry.o :
	$(F90C) -c $(DBLFLAGS) $(F90FLAGS) ../plot/geometry.f90 -o $@

clean :
	rm *.o *.mod
