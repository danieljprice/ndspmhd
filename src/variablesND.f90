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

!-------------------------------------------------------------------
! contains all the global variables used throughout the program
! sorted into modules (alphabetical order)
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! quantities related to the artificial viscosity
!  note alpha and alphain are regarded as particle properties (in module part)
!  similarly the derivative of alpha is in module rates
!-------------------------------------------------------------------
module artvi
 implicit none
 real :: beta,avfact,avdecayconst
 real :: alphamin,alphaumin,alphabmin
end module artvi

!-------------------------------------------------------------------
! boundary related quantities (boundary positions, array storing the 
!   real particle which is mirrored by the ghost)
!-------------------------------------------------------------------

module bound
 use dimen_mhd
 implicit none
 integer, dimension(:), allocatable :: ireal
 real, dimension(ndim) :: xmin, xmax 
 real :: hhmax,pext
end module bound

!-------------------------------------------------------------------
!  debugging quantities
!-------------------------------------------------------------------     

module debug
 implicit none
 integer :: itemp    ! to print debugging info for a specific particle
 logical :: trace
 character(len=6) :: idebug
end module

!-------------------------------------------------------------------
!  curl and divergence of the magnetic field
!-------------------------------------------------------------------

module derivB
 implicit none
 real, dimension(:), allocatable :: divB
 real, dimension(:,:), allocatable :: curlB
end module

!-------------------------------------------------------------------
!  Lorentz force
!-------------------------------------------------------------------

module fmagarray
 implicit none
 real, dimension(:,:), allocatable :: fmag
end module

!-------------------------------------------------------------------
!  correction terms when using a spatially variable smoothing length
!-------------------------------------------------------------------     

module hterms
 implicit none
 integer :: itsdensity
 real, dimension(:), allocatable :: gradh,gradhn,gradsoft,gradgradh,zeta
 real :: rhomin,h_min
end module hterms

!-------------------------------------------------------------------
!  quantities used for link-list neighbour finding
!-------------------------------------------------------------------

module linklist
 use dimen_mhd
 implicit none
 integer, dimension(:), allocatable :: ll,ifirstincell,iamincell,numneigh
 integer, dimension(ndim) :: ncellsx
 integer :: ncells,ncellsloop
 real :: dxcell
end module linklist

!-------------------------------------------------------------------
!  logical unit numbers for input and output files
!-------------------------------------------------------------------

module loguns
 implicit none
 integer :: iprint,ievfile,idatfile,iread,ireadf
 integer :: ifile
 character(len=120) :: rootname       ! name of the run
end module loguns   

!-------------------------------------------------------------------
!  program options
!-------------------------------------------------------------------

module options
 use dimen_mhd
 implicit none
 integer :: iener,icty,iav,ikernav
 integer :: idrag_nature,idrag_structure,ismooth
 integer :: iprterm,idumpghost,ihvar,idust
 integer :: imhd,imagforce,idivbzero   !  (mhd options)
 integer :: iexternal_force,ixsph
 integer :: igravity,ikernel,ikernelalt,iresist
 integer :: maxdensits
 integer :: ivisc,ibiascorrection,iambipolar
 integer, dimension(ndim) :: ibound
 integer, dimension(3) :: iavlim
 real :: damp,dampz,dampr,psidecayfact,tolh,hsoft,etamhd
 real :: Kdrag
 real :: shearvisc,bulkvisc
 real :: gamma_ambipolar,rho_ion
 character(len=12) :: geom
 logical :: usenumdens,use_sqrtdustfrac
end module

!-------------------------------------------------------------------
!  basic particle properties
!-------------------------------------------------------------------

module part
 use dimen_mhd
 implicit none
 integer :: npart,nbpts,ntotal
 integer, parameter :: itypegas = 0
 integer, parameter :: itypebnd = 1
 integer, parameter :: itypedust = 2
 integer, parameter :: itypegas1 = 3
 integer, parameter :: itypegas2 = 4
 integer, parameter :: itypebnd2 = 11
 integer, dimension(:), allocatable :: itype
 real, dimension(:), allocatable    :: pmass,sqrtg
 real, dimension(:,:), allocatable  :: x   
 real, dimension(:), allocatable    :: dens,rho,pr,uu,en,hh,psi,spsound,rhoalt
 real, dimension(:,:), allocatable  :: vel,pmom,sourceterms,alpha
 real, dimension(:,:), allocatable  :: Bfield, Bevol
 real, dimension(ndimB)             :: Bconst
 real, dimension(:,:), allocatable  :: deltav
 real, dimension(:), allocatable    :: dustfrac,dustevol
 real, dimension(:), allocatable    :: del2v
end module part 

!-------------------------------------------------------------------
!  particle properties at the beginning of the time step
!-------------------------------------------------------------------

module part_in
 implicit none
 real, dimension(:), allocatable :: rhoin,prin,hhin,enin,psiin
 real, dimension(:,:), allocatable :: xin,velin,pmomin,alphain
 real, dimension(:,:), allocatable :: Bevolin
 real, dimension(:,:), allocatable :: deltavin
 real, dimension(:),   allocatable :: dustevolin
end module

!-------------------------------------------------------------------
!  rates of change of particle properties (force, derivatives)
!-------------------------------------------------------------------

module rates
 implicit none
 real, dimension(:), allocatable :: drhodt,dudt,dendt,dhdt,dpsidt,poten
 real, dimension(:,:), allocatable :: force,dBevoldt,daldt,gradpsi
 real, dimension(:,:), allocatable :: ddeltavdt
 real, dimension(:), allocatable   :: ddustevoldt
 real :: potengrav
end module rates

!-------------------------------------------------------------------
!  initial particle separation, initial smoothing length (in units of psep)
!-------------------------------------------------------------------

module setup_params
 use dimen_mhd
 use options, only:geom
 implicit none
 character(len=12) :: geomsetup
 real, parameter :: pi = 3.1415926536
!
!--parameters for 2D-MRI simulations
!
 real, parameter :: Rcentre = 1.
 real, parameter :: Omega2 = 1./Rcentre**3
 real, parameter :: Omega0 = Rcentre**(-1.5)
!--for a Keplerian rotation, domegadr = -dlnOmega/dlnr = q = -2A/Omega0 = 1.5 
 real, parameter :: domegadr = 1.5 !--1.5
 real :: psep,hfact, R_grav,xlayer,Alayercs,dwidthlayer
 real, parameter :: omegafixed = 1.0
 
end module setup_params

!-------------------------------------------------------------------
!  time stepping related quantities
!-------------------------------------------------------------------

module timestep
 implicit none
!
!  (max tsteps, number of steps before output, current number of steps,
!   number of steps before using a direct summation)   
!
 integer :: nmax,nout,nsteps,ndirect,nsubsteps_divB,iseedMC
!
! (max time, time before output, dt)
!
 real :: tmax,tout,dt,dt0,time,vsig2max,dtscale
!
! (time step criterion from forces, courant condition)
!
 real :: dtforce, dtcourant, C_force, C_cour, C_rho, dtrho
 real :: dtav, dtdrag, dtvisc
 logical, parameter :: dtfixed = .false.

end module timestep

!-------------------------------------------------------------------
!   version number
!-------------------------------------------------------------------

module versn
 implicit none
 character(len=30) :: version
end module versn

!-------------------------------------------------------------------
! xsph factor
!-------------------------------------------------------------------

module xsph
 implicit none
 real, dimension(:,:), allocatable :: xsphterm
 real :: xsphfac
end module xsph
