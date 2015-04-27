!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

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
 use kernels, only:ianticlump,eps,neps
 implicit none

 if (trace) write(iprint,*) 'Entering subroutine set_default_options'
!
!--set default options
!         
   psep = 0.01
   C_cour = 0.3
   C_force = 0.25
   C_rho = 0.05
   dtforce = 0.
   dtcourant = 0.
   tmax = 1.0
   tout = 0.1
   nmax = 1000000
   nout = -1
   gamma = 5./3.
   iener = 2
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
   tolh = 1.e-3
   idumpghost = 1
   imhd = 0
   imagforce = 2
   idivBzero = 0
   psidecayfact = 0.1
   ianticlump = 0
   eps = 0.8
   neps = 4
   iresist = 0
   etamhd = 0.0
   ixsph = 0
   xsphfac = 0.0
   igravity = 0
   damp = 0.0
   iexternal_force = 0
   ikernel = 0
   ikernelalt = ikernel
   geomsetup = 'cartesian'
   geom = geomsetup
   idim = 0
   usenumdens = .false.
   use_sqrtdustfrac = .false.
   idust = 0
   idrag_nature = 0
   idrag_structure = 0
   Kdrag = 0.
   iambipolar = 0
   rho_ion = 1.e-5
   gamma_ambipolar = 1.
   
 end subroutine set_default_options
