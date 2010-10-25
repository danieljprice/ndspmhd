##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     University of Exeter, UK 2004                                 ##
##-------------------------------------------------------------------##


.KEEP_STATE:

# specify editor to use for file editing
#KNOWN_SYSTEM=no
SHELL = /bin/bash
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

ifndef SETUP1D
#SETUP1D = setup_unifdis1D.f90
#SETUP1D = setup_unifdisND.f90
#SETUP1D = setup_toystar1D_pmass.f90
#SETUP1D = setup_toystar1D.f90
#SETUP1D = setup_shock1DND_mhd2.f90
#SETUP1D = setup_advection1D_mhd.f90
#SETUP1D = setup_wave_x_ND_mhd.f90
#SETUP1D = setup_sphericalblastND_mhd.f90
#SETUP1D = setup_alfvenwaveND_x.f90
SETUP1D = setup_shockND.f90
#SETUP1D = setup_Bxpeak25D_mhd.f90
endif

ifndef SETUP2D
#SETUP2D = setup_orszagtang2D_mhd.f90
#SETUP2D = setup_shock2D_mhd2.f90 
#SETUP2D = setup_sphericalblastND_mhd.f90
#SETUP2D = setup_polytropeND.f90
#SETUP2D = setup_rotor2D_mhd.f90
SETUP2D = setup_unifdisND.f90
#SETUP2D = setup_unifsphND.f90
#SETUP2D = setup_unifdis_cyl.f90
#SETUP2D = setup_toystarND_static.f90
#SETUP2D = setup_shockND.f90
#SETUP2D = setup_shearflow2D.f90
#SETUP2D = setup_alfvenwaveND.f90
#SETUP2D = setup_alfvenwaveND_x.f90
#SETUP2D = setup_Bxpeak25D_mhd.f90
#SETUP2D = setup_kelvinhelmholtz2D.f90
#SETUP2D = setup_blobevap2D.f90
#SETUP2D = setup_cylblast2D.f90
#SETUP2D = setup_jadvect2D.f90
#SETUP2D= setup_nova2D.f90
#SETUP2D= setup_highex2D.f90
endif

ifndef SETUP3D
#SETUP3D = setup_polytropeND.f90
SETUP3D = setup_unifdisND.f90
#SETUP3D = setup_unifsphND.f90
#SETUP3D = setup_sphericalblastND_mhd.f90
#SETUP3D = setup_wave_x_ND_mhd.f90
#SETUP3D = setup_shockND.f90
#SETUP3D = setup_alfvenwaveND_x.f90
#SETUP3D = setup_densityprofileND.f90 
#SETUP3D = setup_orszagtang2D_mhd.f90
endif

######################################################################
########### setup file for modify dump utility #######################

MODDUMPSRC1D = moddump_tstar2D.f90
MODDUMPSRC2D = moddump_tstar2D.f90
MODDUMPSRC3D = moddump_tstar2D.f90
######################################################################

CONS2PRIM = conservative2primitive.f90

# these files are common to all dimensions and coordinate systems

SOURCES= ndspmhd.f90 eos.f90 \
         allocateND.f90 kernelND.f90 readwrite_dumps.f90 \
         get_neighbour_lists.f90 smooth.f90 \
         readwrite_infile.f90 set_uniform_distributionND.f90 \
         boundaryND.f90 check_setup.f90 get_curl.f90 ${CONS2PRIM} \
	 copy_particle.f90 defaults.f90 derivs.f90 \
         direct_sum_poisson3D_soft.f90 \
         density_sums.f90 evolve.f90 external_forces.f90 \
         ghostND_mhd.f90 initialiseND_mhd.f90 \
	 iterate_density.f90 linkND.f90 logun.f90 outputND_mhd.f90\
	 random.f90 ratesND_mhd.f90 riemannsolver.f90 \
	 set_fixedbound.f90 utils.f90 write_headerND_mhd.f90 \
	 get_divB.f90 divBcorrect.f90 \
	 check_neighbourlist.f90 direct_sum_poisson2D_vec.f90

# these are the ordinary cartesian versions
SOURCES2 = stepND_mhd.f90 evwrite_mhd.f90 
	 
OBJECTS1D =  $(SOURCES1D:.f90=.o) $(SOURCES:.f90=.o) $(SOURCES2:.f90=.o) \
	    $(SETUP1D:.f90=.o)
OBJECTS2D =  $(SOURCES2D:.f90=.o) $(SOURCES:.f90=.o)  $(SOURCES2:.f90=.o) \
	    $(SETUP2D:.f90=.o)
OBJECTS3D =  $(SOURCES3D:.f90=.o) $(SOURCES:.f90=.o) $(SOURCES2:.f90=.o)  \
	    $(SETUP3D:.f90=.o)

crap:
	echo specify number of dimensions, e.g. make 1D

1D: checksystem $(OBJECTS1D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../1DSPMHD $(OBJECTS1D)

2D : checksystem $(OBJECTS2D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../2DSPMHD $(OBJECTS2D)

3D : checksystem $(OBJECTS3D)
	$(F90C) $(F90FLAGS) $(LDFLAGS) -o ../3DSPMHD $(OBJECTS3D)

## modify dump utility

OBJMODDUMP = variablesND.o allocateND.o \
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
	@${MAKE} compilers
	@$(MAKE) err;
   endif

compilers:
	@echo "I suggest one of the following, based on detected Fortran compilers..."; echo;
	@if type -p ifort > /dev/null; then echo "make SYSTEM=ifort"; fi;
	@if type -p pathf90 > /dev/null; then echo "make SYSTEM=pathf90"; fi;
	@if type -p pgf90 > /dev/null; then echo "make SYSTEM=pgf90"; fi;
	@if type -p xlf90_r > /dev/null; then echo "make SYSTEM=ukaff1a [uses xlf90_r]"; fi;
	@if type -p gfortran > /dev/null; then echo "make SYSTEM=gfortran"; fi;
	@if type -p g95 > /dev/null; then echo "make SYSTEM=g95"; fi;
	@echo "(end of possible selections)"; echo;

err:
	$(error aborting);

## other crap

tar:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90

targz:
	tar cf NDSPMHD.tar Makefile dimen*.f90 $(SOURCES) setup*.f90
	gzip NDSPMHD.tar

edit1D: checkeditor
	$(EDITOR) $(SETUP1D) &
edit2D: checkeditor
	$(EDITOR) $(SETUP2D) &
edit3D: checkeditor
	$(EDITOR) $(SETUP3D) &

checkeditor:
   ifndef EDITOR
	@echo
	@echo 'ERROR: need EDITOR variable set for "make edit" (e.g. export EDITOR=emacs)'
	@echo
	@$(MAKE) err;
   endif

clean:
	rm *.o *.mod
