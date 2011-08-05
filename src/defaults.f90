!!-----------------------------------------------------------------
!! Sets default values of input parameters
!! these are overwritten by reading from the input file or 
!! by setting them in the setup routine
!!-----------------------------------------------------------------

subroutine set_default_options
 use dimen_mhd
 use debug
 use loguns
 
 use artvi
 use eos
 use options
 use setup_params
 use timestep
 use xsph
 use anticlumping
 implicit none

 if (trace) write(iprint,*) 'Entering subroutine set_default_options'
!
!--set default options
!         
   psep = 0.01
   C_cour = 0.3
   C_force = 0.25
   dtforce = 0.
   dtcourant = 0.
   tmax = 1.0
   tout = 0.01
   nmax = 1000000
   nout = -1
   gamma = 5./3.
   iener = 3
   polyk = 1.0
   icty = 0
   ndirect = nmax
   maxdensits = 250
   iprterm = 0
   iav = 2
   alphamin = 0.1
   alphaumin = 0.0
   alphaBmin = 1.0
   beta = 2.0
   iavlim(1) = 2
   iavlim(2) = 1
   iavlim(3) = 0
   avdecayconst = 0.1
   ikernav = 3
   ihvar = 2
   hfact = 1.2
   idumpghost = 1
   imhd = 1
   imagforce = 2
   idivBzero = 0
   psidecayfact = 0.1
   ianticlump = 0
   eps = 0.8
   neps = 4
   ixsph = 0
   xsphfac = 0.0
   igravity = 0
   damp = 0.0
   iexternal_force = 0
   ikernel = 0
   igeom = 1
   igeomsetup = 1
   
 end subroutine set_default_options
