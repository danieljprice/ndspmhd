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
 use polyconst
 use setup_params
 use timestep
 use xsph
 use anticlumping

 if (trace) write(iprint,*) 'Entering subroutine set_default_options'
!
!--set default options
!         
   psep = 0.01
   C_cour = 0.8
   C_force = 0.25
   tmax = 1.0
   tout = 0.01
   nmax = 1000000
   nout = -1
   gamma = 5./3.
   iener = 3
   udiss_frac = 1.0
   Bdiss_frac = 1.0
   polyk = 1.0
   icty = 0
   ndirect = nmax
   ialtform = 0
   iav = 1
   alphamin = 0.1
   beta = 2.0
   iavlim = 1
   avconst = 0.1
   ikernav = 3
   ihvar = 3
   hfact = 1.3
   idumpghost = 1
   imhd = 1
   imagforce = 2
   idivBzero = 0
   ianticlump = 1
   eps = 0.8
   neps = 5
   ixsph = 0
   xsphfac = 0.0
   igravity = 0
   damp = 0.0
   iexternal_force = 0
   ikernel = 0
   
 end subroutine set_default_options
