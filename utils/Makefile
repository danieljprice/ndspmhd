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
FC =  gfortran
FFLAGS =  -O2 -Wall -fdefault-real-8 -fdefault-double-8 #-fbounds-check
#LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot95 -lpgplot
LDFLAGS =  -L/usr/X11/lib -lX11 -L${PGPLOT_DIR} -lpgplot -lpng
DBLFLAGS = -fdefault-real-8
ifndef SPLASH_DIR
SPLASH_DIR=$(HOME)/splash/src
endif

.SUFFIXES:
.SUFFIXES:.f .f90 .o .mod
VPATH=../src $(SPLASH_DIR)
CC = $(FC)

# define the implicit rule to make a .o file from a .f90 file
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# define the implicit rule for compiling directly from a single .f90 file
#%.f90 :
#	$(FC) $(FFLAGS) $(LDFLAGS) $< -o $@

OBJECTS=kernel_cspline_utils.o kernel_erbs_utils.o kernelND.o

OBJSPHDERIV = kernelND.o sphderiv.o

OBJMODDUMP = variablesND.o allocateND.o prompting.o \
             readwrite_dumps.o modifydump.o

OBJPLOT = legend.o setpage.o

all: kernelplot kernelplot3D kernelcalc3D getkernellist

kernelplot: $(OBJECTS) $(OBJPLOT) kernelplot.o kernelstability1D.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(OBJPLOT) kernelplot.o kernelstability1D.o $(LDFLAGS)

kernelplot3D: $(OBJECTS) $(OBJPLOT) kernel_sums.o kernelplot3D.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(OBJPLOT) kernel_sums.o kernelplot3D.o $(LDFLAGS)

kernelcalc3D: $(OBJECTS) kernel_sums.o kernelcalc3D.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) kernel_sums.o kernelcalc3D.o

getkernellist: $(OBJECTS) getkernellist.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) getkernellist.o $(LDFLAGS)

SRCSORT=utils_sort.f90 sort.f90
OBJSORT=${SRCSORT:.f90=.o}
sort: $(OBJSORT)
	$(FC) $(FFLAGS) -o $@ $(OBJSORT)

SRCBREED=random.f90 breed.f90
OBJBREED=${SRCBREED:.f90=.o}
breed: $(OBJBREED)
	$(FC) $(FFLAGS) -o $@ $(OBJBREED)

kernelplot.o : kernelND.o

sametime : sametime.o
	$(FC) $(FFLAGS) -o $@ sametime.o $(LDFLAGS)

plotfreq : plotfreq.o danpgtile.o
	$(FC) $(FFLAGS) -o $@ plotfreq.o danpgtile.o $(LDFLAGS)

plotfreq.o : plotfreq.f90
	$(FC) -c $(FFLAGS) plotfreq.f90 -o $@

sphderiv : $(OBJSPHDERIV)
	$(FC) $(FFLAGS) -o $@ $(OBJSPHDERIV) $(LDFLAGS)
        
moddump1D: dimen1D.o $(OBJMODDUMP)
	$(FC) $(FFLAGS) -o $@ dimen1D.o $(OBJMODDUMP)

moddump2D: dimen2D.o $(OBJMODDUMP)
	$(FC) $(FFLAGS) -o $@ dimen2D.o $(OBJMODDUMP)

moddump3D: dimen3D.o $(OBJMODDUMP)
	$(FC) $(FFLAGS) -o $@ dimen3D.o $(OBJMODDUMP)

dimen1D.o :
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/dimen1D.f90 -o $@

dimen2D.o :
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/dimen2D.f90 -o $@

dimen3D.o :
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/dimen3D.f90 -o $@

variablesND.o :
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/variablesND.f90 -o $@

kernelND.o :
	$(FC) -c $(FFLAGS) ../src/kernelND.f90 -o $@

allocateND.o : variablesND.o
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/allocateND.f90 -o $@

readwrite_dumps.o : variablesND.o geometry.o
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../src/readwrite_dumps.f90 -o $@

modifydump.o : modifydump.f90
	$(FC) -c $(DBLFLAGS) $(FFLAGS) modifydump.f90 -o $@

geometry.o :
	$(FC) -c $(DBLFLAGS) $(FFLAGS) ../plot/geometry.f90 -o $@

clean :
	rm *.o *.mod
