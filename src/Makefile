##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     University of Exeter, UK 2004                                 ##
##-------------------------------------------------------------------##


.KEEP_STATE:

# specify editor to use for file editing
EDITOR = nedit
KNOWN_SYSTEM=no
SHELL = /bin/sh

DOUBLEPRECISION= yes

ifeq ($(SYSTEM), ukaff1a)
    F90C= xlf_r -qnosave
    F90FLAGS= -O3 -q64 -qextname -qsmp=noauto
    DBLFLAG= -qrealsize=8
#    DEBUGFLAG= -qflttrap=enable:invalid:zerodivide -g -C -qsigtrap -qfloat=nans 
    DEBUGFLAG= -C
    KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), sunf95)
    F90C= sunf95
    F90FLAGS= -fast
    DBLFLAG= -xtypemap=real:64,double:64 -ftrap=common
    DEBUGFLAG= -g -C -ftrap=common
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -xfilebyteorder=big16:%all ##:45
    ENDIANFLAGLITTLE= -xfilebyteorder=little16:%all
endif

ifeq ($(SYSTEM), g95)
    F90C= g95
    F90FLAGS= -O3 -Wall -Wno=112,140 -ffast-math -funroll-loops -finline-functions -falign-loops=16
    DBLFLAG= -r8
    DEBUGFLAG= -g #-ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
# or use setenv G95_ENDIAN='BIG' or 'LITTLE' at runtime
endif

ifeq ($(SYSTEM), mymac)
    F90C= g95
    F90FLAGS= -O3 -g -Wall -Wno=112,140 -ffast-math -funroll-loops -finline-functions -falign-loops=16
    DBLFLAG= -r8
    DEBUGFLAG= -g -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
# or use setenv G95_ENDIAN='BIG' or 'LITTLE' at runtime
endif

ifeq ($(SYSTEM), gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   F90C= gfortran
   F90FLAGS= -O3 -Wall -funroll-loops -g #-ffast-math #-mtune=native -ftree-vectorize -ftree-vectorizer-verbose=3
   DBLFLAG= -fdefault-real-8
   SYSTEMFILE= system_unix.f90
   DEBUGFLAG= -g -frange-check
   OMPFLAG= -fopenmp
   ENDIANFLAGBIG= -fconvert=big-endian
   ENDIANFLAGLITTLE= -fconvert=little-endian
   KNOWN_SYSTEM=yes
endif

ifeq ($(SYSTEM), myg95)
    F90C= myg95
    F90FLAGS= -O3 -Wall -Wno=112,140 -Wextra -ffast-math 
    DBLFLAG= -r8
    DEBUGFLAG= -ftrace=full -fbounds-check
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -fendian='BIG'
    ENDIANFLAGLITTLE= -fendian='LITTLE'
endif

ifeq ($(SYSTEM), ifort)
    F90C= ifort
    F90FLAGS= -O3 -Vaxlib -w90 -w95 -cm
    DBLFLAG= -r8
    DEBUGFLAG= -C
    KNOWN_SYSTEM=yes
    ENDIANFLAGBIG= -convert big_endian
    ENDIANFLAGLITTLE= -convert little_endian
# or use setenv F_UFMTENDIAN=big:45 at runtime (e.g. for unit 45 only)
endif


# Set other optional flags depending on settings

ifeq ($(DEBUG), yes)
    F90FLAGS += ${DEBUGFLAG}
endif

ifeq ($(ENDIAN), BIG)
    F90FLAGS += ${ENDIANFLAGBIG}
endif

ifeq ($(ENDIAN), LITTLE)
    F90FLAGS += ${ENDIANFLAGLITTLE}
endif

ifeq ($(DOUBLEPRECISION), yes)
    F90FLAGS += ${DBLFLAG}
endif

# Fortran flags same as F90
FC = $(F90C)
FFLAGS = $(F90FLAGS)

# define the implicit rule to make a .o file from a .f90 file

%.o : %.f90
	$(F90C) -c $(F90FLAGS) $(FPPFLAGS) $< -o $@


# use the appropriate dimensions module and dimension specific routines

SOURCES1D = dimen1D.f90 variablesND.f90 direct_sum_poisson1D.f90
SOURCES2D = dimen2D.f90 variablesND.f90 direct_sum_poisson2D.f90
SOURCES3D = dimen3D.f90 variablesND.f90 direct_sum_poisson3D.f90 

######################################################################
########### choose the setup file in each dimension here #############

#SETUP1D = setup_unifdis1D.f90
#SETUP1D = setup_unifdisND.f90
#SETUP1D = setup_toystar1D_pmass.f90
#SETUP1D = setup_toystar1D.f90
#SETUP1D = setup_shock1DND_mhd2.f90
#SETUP1D = setup_advection1D_mhd.f90
#SETUP1D = setup_wave_x_ND_mhd.f90
#SETUP1D = setup_turbulenceND_mhd.f90 set_powerspec1D.f90 
#SETUP1D = setup_sphericalblastND_mhd.f90
#SETUP1D = setup_alfvenwaveND_x.f90
SETUP1D = setup_shockND.f90
#SETUP1D = setup_acousticshockND.f90
#SETUP1D = setup_Bxpeak25D_mhd.f90
#SETUP1D = setup_clareshk1D.f90
#SETUP1D = setup_resist1D.f90
#SETUP1D = setup_anyrho1Dran.f90

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
#SETUP2D = setup_shockND.f90
#SETUP2D = setup_fromfile.f90
#SETUP2D = setup_shearflow2D.f90
#SETUP2D = setup_rings2D.f90
#SETUP2D = setup_alfvenwaveND.f90
#SETUP2D = setup_alfvenwaveND_x.f90
#SETUP2D = setup_Bxpeak25D_mhd.f90
#SETUP2D = setup_disc2D_2.f90
#SETUP2D = setup_logodivB.f90
#SETUP2D = setup_mri2D.f90
#SETUP2D = setup_torus2D.f90
SETUP2D = setup_kelvinhelmholtz2D.f90
#SETUP2D = setup_blobevap2D.f90
#SETUP2D = setup_cylblast2D.f90
#SETUP2D = setup_jadvect2D.f90
#SETUP2D = setup_neutronstarkh.f90
#SETUP2D = setup_sunblast2D.f90
#SETUP2D = setup_recontest3D.f90
#SETUP2D = setup_convection2D.f90
#SETUP2D = setup_anyrho2Dmin.f90
#SETUP2D = setup_posfromfile.f90
#SETUP2D = setup_rayleightaylor2D.f90
#SETUP2D = setup_blobweirdshape.f90
#SETUP2D= setup_taylorcouette2D.f90
#SETUP2D= setup_nova2D.f90
#SETUP2D= setup_highex2D.f90

#SETUP3D = setup_polytropeND.f90
#SETUP3D = setup_unifdisND.f90 #cp_distribute.f
#SETUP3D = setup_unifsphND.f90
#SETUP3D = setup_sphericalblastND_mhd.f90
#SETUP3D = setup_wave_x_ND_mhd.f90
SETUP3D = setup_shockND.f90
#SETUP3D = setup_alfvenwaveND_x.f90
#SETUP3D = setup_disc3D.f90
#SETUP3D = setup_torus3D.f90
#SETUP3D = setup_fromfile.f90
#SETUP3D = setup_densityprofileND.f90 
#SETUP3D = setup_recontest3D.f90
#SETUP3D = setup_forcefree3D.f90
#SETUP3D = setup_blastisoND.f90
#SETUP3D = setup_resistND.f90

######################################################################
########### setup file for modify dump utility #######################

MODDUMPSRC1D = moddump_tstar2D.f90
MODDUMPSRC2D = moddump_tstar2D.f90
MODDUMPSRC3D = moddump_tstar2D.f90
######################################################################

CONS2PRIM = conservative2primitive.f90

# these files are common to all dimensions and coordinate systems

SOURCES= sphND_mhd.f90 geometry.f90 convert_gr.f90 eos.f90 \
         allocateND.f90 kernelND.f90 readwrite_dumps.f90 \
         get_neighbour_lists.f90 \
         readwrite_infile.f90 set_uniform_distributionND.f90 \
         boundaryND.f90 check_setup.f90 get_curl.f90 smooth.f90 ${CONS2PRIM} \
	 copy_particle.f90 particlesplitting.f90 \
	 defaults.f90 derivs.f90 \
         direct_sum_poisson3D_soft.f90 \
         grutils_gr.f90 \
         density_sums.f90 densityiterate.f90 \
	 exact_densityprofiles.f90 errors.f90 \
         evolve.f90 external_forces.f90 \
         ghostND_mhd.f90 initialiseND_mhd.f90 \
	 iterate_density.f90 \
	 linkND.f90 logun.f90 \
         outputND_mhd.f90\
	 random.f90 ratesND_mhd.f90 \
	 riemannsolver.f90 \
	 set_fixedbound.f90 \
	 sort_particles.f90 substep_divB.f90 utils.f90 \
	 write_headerND_mhd.f90 \
	 get_divB.f90 divBcorrect.f90 get_divBgradpsi.f90 \
	 check_neighbourlist.f90 direct_sum_poisson2D_vec.f90

# these are the ordinary cartesian versions
SOURCES2 = stepND_leapfrog_mhd.f90 evwrite_mhd.f90 

# these files are for the sr/gr version of the code
SOURCESGR = stepND_mhd_gr.f90 evwrite_mhd_gr.f90
CONS2PRIMGR = conservative2primitive_sr.f90
	 
OBJECTS1D =  $(SOURCES1D:.f90=.o) $(SOURCES:.f90=.o) $(SOURCES2:.f90=.o) \
	    $(SETUP1D:.f90=.o)
OBJECTS2D =  $(SOURCES2D:.f90=.o) $(SOURCES:.f90=.o)  $(SOURCES2:.f90=.o) \
	    $(SETUP2D:.f90=.o)
OBJECTS3D =  $(SOURCES3D:.f90=.o) $(SOURCES:.f90=.o) $(SOURCES2:.f90=.o)  \
	    $(SETUP3D:.f90=.o)
            
SOURCES1GR = $(SOURCES:${CONS2PRIM}=${CONS2PRIMGR})

OBJECTSGR1D = $(SOURCES1D:.f90=.o) $(SOURCES1GR:.f90=.o) $(SOURCESGR:.f90=.o)\
	    $(SETUP1D:.f90=.o)
OBJECTSGR2D = $(SOURCES2D:.f90=.o) $(SOURCES1GR:.f90=.o) $(SOURCESGR:.f90=.o) \
	    $(SETUP2D:.f90=.o)
OBJECTSGR3D = $(SOURCES3D:.f90=.o) $(SOURCES1GR:.f90=.o) $(SOURCESGR:.f90=.o) \
	    $(SETUP3D:.f90=.o)

crap:
	echo specify number of dimensions, e.g. make 1D

1D: checksystem $(OBJECTS1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DSPMHD $(OBJECTS1D)

2D : checksystem $(OBJECTS2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DSPMHD $(OBJECTS2D)

3D : checksystem $(OBJECTS3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DSPMHD $(OBJECTS3D)

1DGR : checksystem $(OBJECTSGR1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DGRSPMHD $(OBJECTSGR1D)

2DGR : checksystem $(OBJECTSGR2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DGRSPMHD $(OBJECTSGR2D)

3DGR : checksystem $(OBJECTSGR3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DGRSPMHD $(OBJECTSGR3D)

## files from elsewhere

prompting.o: ../plot/prompting.f90
	$(F90C) -c ../plot/prompting.f90 -o $@

geometry.o: ../plot/geometry.f90
	$(F90C) -c $(F90FLAGS) ../plot/geometry.f90 -o $@

exact_densityprofiles.o: ../plot/exact_densityprofiles.f90
	$(F90C) -c $(F90FLAGS) ../plot/exact_densityprofiles.f90 -o $@

evolve_montecarlo.o: setup_densityprofileND.o

## modify dump utility

OBJMODDUMP = variablesND.o allocateND.o geometry.o \
             readwrite_dumps.o prompting.o modifydump.o
OBJMODDUMP1D = $(OBJMODDUMP) $(MODDUMPSRC1D:.f90=.o)
OBJMODDUMP2D = $(OBJMODDUMP) $(MODDUMPSRC2D:.f90=.o)
OBJMODDUMP3D = $(OBJMODDUMP) $(MODDUMPSRC3D:.f90=.o)
             

moddump1D: dimen1D.o $(OBJMODDUMP1D)
	$(F90C) $(F90FLAGS) -o ../$@ dimen1D.o $(OBJMODDUMP1D)

moddump2D: dimen2D.o $(OBJMODDUMP2D)
	$(F90C) $(F90FLAGS) -o ../$@ dimen2D.o $(OBJMODDUMP2D)

moddump3D: dimen3D.o $(OBJMODDUMP3D)
	$(F90C) $(F90FLAGS) -o ../$@ dimen3D.o $(OBJMODDUMP3D)

checksystem:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling ndspmhd for $(SYSTEM) system..........."
	@echo ""
        ifeq ($(DEBUG), yes)
	     @echo "Debugging flags are ON"
        endif
        ifeq ($(DOUBLEPRECISION), yes)
	     @echo "Flags set for DOUBLE PRECISION"
        else
	     @echo "Flags set for SINGLE PRECISION"
        endif
        ifeq ($(ENDIAN), BIG)
	     @echo "Flags set for conversion to BIG endian"
        endif
        ifeq ($(ENDIAN), LITTLE)
	     @echo "Flags set for conversion to LITTLE endian"
        endif
   else
	@echo ""
	@echo "make: WARNING: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=>set the environment variable SYSTEM to one listed "
	@echo "  in the Makefile and try again"
	@echo ""
	quit
   endif

## other crap

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
