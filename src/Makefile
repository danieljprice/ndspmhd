##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2003	 	     ##
##-------------------------------------------------------------------##


.KEEP_STATE:

# specify editor to use for file editing
EDITOR = nedit

#F90C = ifc
#F90FLAGS =  -O -C -WB -r8
#LDFLAGS = -Vaxlib
F90C =  g95
F90FLAGS = -O -r8 -Wall##-fbounds-check -Wall
LDFLAGS =  -L/usr/X11R6/lib -lX11 -L/sw/lib -lpng -laquaterm -lcc_dynamic -Wl,-framework -Wl,Foundation -L/sw/lib/pgplot95 -lpgplot

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@


# use the appropriate dimensions module

SOURCES1D = dimen1D.f90 variablesND.f90 direct_sum_poisson1D.f90
SOURCES2D = dimen2D.f90 variablesND.f90 direct_sum_poisson2D.f90
SOURCES3D = dimen3D.f90 variablesND.f90 direct_sum_poisson3D.f90

######################################################################
########### choose the setup file in each dimension here #############

#SETUP1D = setup_unifdis1D.f90
#SETUP1D = setup_unifdisND.f90
#SETUP1D = setup_toystar1D_pmass.f90
#SETUP1D = setup_toystar1D.f90
SETUP1D = setup_shock1DND_mhd2.f90
#SETUP1D = setup_advection1D_mhd.f90
#SETUP1D = setup_wave_x_ND_mhd.f90
#SETUP1D = setup_turbulenceND_mhd.f90 set_powerspec1D.f90 
#SETUP1D = setup_sphericalblastND_mhd.f90
#SETUP1D = setup_alfvenwaveND_x.f90
#SETUP1D = setup_shockND.f90
#SETUP1D = setup_acousticshockND.f90
#SETUP1D = setup_Bxpeak25D_mhd.f90

#SETUP2D = setup_orszagtang2D_mhd.f90
#SETUP2D = setup_shock2D_mhd2.f90 
#SETUP2D = setup_sphericalblastND_mhd.f90
#SETUP2D = setup_polytropeND.f90
#SETUP2D = setup_rotor2D_mhd.f90
#SETUP2D = setup_unifdisND.f90
#SETUP2D = setup_unifsphND.f90
#SETUP2D = setup_unifdis_cyl.f90
#SETUP2D = setup_disc2D.f90
#SETUP2D = setup_toystarND_static.f90 
SETUP2D = setup_shockND.f90
#SETUP2D = setup_fromfile.f90
#SETUP2D = setup_shearflow2D.f90
#SETUP2D = setup_rings2D.f90
#SETUP2D = setup_alfvenwaveND.f90
#SETUP2D = setup_alfvenwaveND_x.f90
#SETUP2D = setup_Bxpeak25D_mhd.f90
#SETUP2D = setup_disc2D_2.f90
#SETUP2D = setup_logodivB.f90

#SETUP3D = setup_unifdisND.f90
#SETUP3D = setup_sphericalblastND_mhd.f90
#SETUP3D = setup_wave_x_ND_mhd.f90
#SETUP3D = setup_shockND.f90
#SETUP3D = setup_alfvenwaveND_x.f90
SETUP3D = setup_disc3D.f90
######################################################################

# these files are common to all dimensions and coordinate systems

SOURCES = sphND_mhd.f90 geometry.f90 convert_gr.f90 \
         allocateND.f90 readwrite_dumps.f90 \
         readwrite_infile.f90 set_uniform_distributionND.f90 \
         boundaryND_1D.f90 check_setup.f90 \
	 copy_particle.f90 \
	 defaults.f90 \
         get_neighbour_lists.f90 \
         density_sums.f90 \
	 eos.f90 evolve.f90 \
	 external_forces.f90 \
         ghostND_mhd.f90 initialiseND_mhd.f90 \
	 interpolate_kernel.f90\
	 iterate_density.f90 kernelND.f90 \
	 linkND.f90 logun.f90 \
         outputND_mhd.f90 \
	 random.f90 ratesND_mhd.f90 \
	 riemannsolver.f90 \
	 set_fixedbound.f90 \
	 sobolsequence.f90 substep_divB.f90 utils.f90 \
	 write_headerND_mhd.f90 \
	 get_divB.f90 divBcorrect.f90 get_curl.f90 get_divBgradpsi.f90 \
	 check_neighbourlist.f90 direct_sum_poisson2D_vec.f90

# these are the ordinary cartesian versions
SOURCES2 = conservative2primitive.f90 stepND_mhd.f90 evwrite_mhd.f90 \

# these files are for the gr version of the code
SOURCESGR = grutils_gr.f90 conservative2primitive_gr.f90 \
	    stepND_mhd_gr.f90 evwrite_mhd_gr.f90
	 
OBJECTS1D =  $(SOURCES1D:.f90=.o) $(SOURCES2:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP1D:.f90=.o)
OBJECTS2D =  $(SOURCES2D:.f90=.o) $(SOURCES2:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP2D:.f90=.o)
OBJECTS3D =  $(SOURCES3D:.f90=.o) $(SOURCES2:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP3D:.f90=.o)

OBJECTSGR1D = $(SOURCES1D:.f90=.o) $(SOURCESGR:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP1D:.f90=.o)

OBJECTSGR2D = $(SOURCES2D:.f90=.o) $(SOURCESGR:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP2D:.f90=.o)

OBJECTSGR3D = $(SOURCES3D:.f90=.o) $(SOURCESGR:.f90=.o) $(SOURCES:.f90=.o) \
	    $(SETUP3D:.f90=.o)

crap:
	echo specify number of dimensions, e.g. make 1D

1D: $(OBJECTS1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DSPMHD $(OBJECTS1D)

2D : $(OBJECTS2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DSPMHD $(OBJECTS2D)

3D : $(OBJECTS3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DSPMHD $(OBJECTS3D)

1DGR : 	$(OBJECTSGR1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DGRSPMHD $(OBJECTSGR1D)

2DGR : 	$(OBJECTSGR2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DGRSPMHD $(OBJECTSGR2D)

3DGR : 	$(OBJECTSGR3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DGRSPMHD $(OBJECTSGR3D)

## files from elsewhere

geometry.o:
	$(F90C) -c $(F90FLAGS) ../plot/geometry.f90 -o $@

tar:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90

targz:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90
	gzip NDSPMHD.tar

edit1D:
	$(EDITOR) $(SETUP1D) &
edit2D:
	$(EDITOR) $(SETUP2D) &
edit3D:
	$(EDITOR) $(SETUP3D) &
clean:
	rm *.o *.mod
