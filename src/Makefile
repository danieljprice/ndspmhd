##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2003	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

FORTRAN = f90
FC = f90
EDITOR = nedit

# -C gives array out-of-bounds checking
# -u forces implicit none
# -g for debugging with dbx etc
# -pg for profiling
# -xtypemap=real:64,double:64,integer:32 for double precision
# -ftrap = %all traps floating point exceptions

#F90FLAGS = -O4 
# use the following instead for double precision
F90FLAGS = -C -O4 -xtypemap=real:64,double:64,integer:32
# PGPLOT libraries (note linking with the -lF77 libraries causes an 
# error in the flush commands, just means output is buffered).
# NOTE: cannot use the double precision flags when compiling with PGPLOT
#LDFLAGS = -lpgplot -lX11 -lF77
LDFLAGS =

# use the appropriate dimensions module

SOURCES1D = dimen1D.f90 
SOURCES2D = dimen2D.f90
SOURCES25D = dimen25D.f90
SOURCES3D = dimen3D.f90

######################################################################
########### choose the setup file in each dimension here #############

#SETUP1D = setup_unifdis1D.f90
#SETUP1D = setup_unifdisND.f90
#SETUP1D = setup_toystar1D_pmass.f90
SETUP1D = setup_toystar1D.f90
#SETUP1D = setup_shock1DND_mhd2.f90
#SETUP1D = setup_advection1D_mhd.f90
#SETUP1D = setup_wave_x_ND_mhd.f90
#SETUP1D = setup_turbulenceND_mhd.f90	 
#SETUP1D = setup_sphericalblastND_mhd.f90

SETUP2D = setup_orszagtang2D_mhd.f90
#SETUP2D = setup_shock2D_mhd2.f90 
#SETUP2D = setup_sphericalblastND_mhd.f90
#SETUP2D = setup_polytropeND.f90
#SETUP2D = setup_rotor2D_mhd.f90
#SETUP2D = setup_unifdisND.f90

SETUP25D = setup_alfvenwaveND.f90

#SETUP3D = setup_unifdisND.f90
#SETUP3D = setup_sphericalblastND_mhd.f90
SETUP3D = setup_shockND.f90

######################################################################

# these files are common to all dimensions

SOURCES = variablesND.f90 sphND_mhd.f90   \
         allocateND.f90 boundaryND_1D.f90 \
	 conservative2primitive.f90 densityND.f90 \
	 densityND_partial.f90 \
	 direct_sum_grav.f90 \
	 eos.f90 evwrite_mhd.f90 get_neighbour_listND.f90 \
	 get_neighbour_listND_partial.f90 \
         ghostND_mhd.f90 initialiseND_mhd.f90 \
	 interpolate_kernel.f90 \
	 iterate_density.f90 kernelND.f90 \
	 linkND.f90 logun.f90 outputND_mhd.f90 \
	 random.f90 \
	 ratesND_mhd.f90 \
	 read_infile_mhd.f90 \
	 set_uniform_cartesianND.f90 set_uniform_sphericalND.f90 \
	 set_velfield1D.f90 \
	 smooth.f90 stepND_leapfrog_mhd.f90 \
	 write_infile_mhd.f90 write_headerND_mhd.f90 \
#	 check_neighbourlist.f90
#	 ratesND_mhd_split.f90 rates_artvisc.f90 rates_energy.f90 rates_mhdterms.f90 	 \
	 
OBJECTS1D = $(SOURCES1D:.f90=.o) $(SOURCES:.f90=.o) $(SETUP1D:.f90=.o)
OBJECTS2D = $(SOURCES2D:.f90=.o) $(SOURCES:.f90=.o) $(SETUP2D:.f90=.o)
OBJECTS25D = $(SOURCES25D:.f90=.o) $(SOURCES:.f90=.o) $(SETUP25D:.f90=.o)
OBJECTS3D = $(SOURCES3D:.f90=.o) $(SOURCES:.f90=.o) $(SETUP3D:.f90=.o)

1D: $(OBJECTS1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DSPMHD $(OBJECTS1D)
	
2D : $(OBJECTS2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DSPMHD $(OBJECTS2D)

25D : $(OBJECTS25D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../25DSPMHD $(OBJECTS25D)

3D : $(OBJECTS3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DSPMHD $(OBJECTS3D)
		
tar:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90

targz:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90
	gzip NDSPMHD.tar

edit1D:
	$(EDITOR) $(SETUP1D) &
edit2D:
	$(EDITOR) $(SETUP2D) &
edit25D:
	$(EDITOR) $(SETUP25D) &
edit3D:
	$(EDITOR) $(SETUP3D) &
clean:
	rm *.o *.mod
