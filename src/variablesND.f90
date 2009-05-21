!-------------------------------------------------------------------
! contains all the global variables used throughout the program
! sorted into modules (alphabetical order)
!-------------------------------------------------------------------

!-------------------------------------------------------------------
! quantities related to the artificial viscosity
!  note alpha and alphain are regarded as particle properties (in module part)
!  similarly the derivative of alpha is in module rates
!-------------------------------------------------------------------
MODULE artvi
 IMPLICIT NONE
 REAL :: beta,avfact,avdecayconst
 REAL :: alphamin,alphaumin,alphaBmin
END MODULE artvi

!-------------------------------------------------------------------
! boundary related quantities (boundary positions, array storing the 
!   real particle which is mirrored by the ghost)
!-------------------------------------------------------------------

MODULE bound
 USE dimen_mhd
 IMPLICIT NONE
 INTEGER, DIMENSION(:), ALLOCATABLE :: ireal
 REAL, DIMENSION(ndim) :: xmin, xmax 
 REAL :: hhmax,pext
END MODULE bound

!-------------------------------------------------------------------
!  debugging quantities
!-------------------------------------------------------------------     

MODULE debug
 IMPLICIT NONE
 INTEGER :: itemp    ! to print debugging info for a specific particle
 LOGICAL :: trace
 CHARACTER(LEN=6) :: idebug
END MODULE

!-------------------------------------------------------------------
!  curl and divergence of the magnetic field
!-------------------------------------------------------------------

MODULE derivB
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: divB
 REAL, DIMENSION(:,:), ALLOCATABLE :: curlB
END MODULE

!-------------------------------------------------------------------
!  Lorentz force
!-------------------------------------------------------------------

MODULE fmagarray
 IMPLICIT NONE
 REAL, DIMENSION(:,:), ALLOCATABLE :: fmag
END MODULE

!-------------------------------------------------------------------
!  correction terms when using a spatially variable smoothing length
!-------------------------------------------------------------------     

MODULE hterms
 IMPLICIT NONE
 INTEGER :: itsdensity
 REAL, DIMENSION(:), ALLOCATABLE :: gradh,gradhn,gradsoft,gradgradh,zeta
 REAL :: rhomin
END MODULE hterms

!-------------------------------------------------------------------
!  correction terms to make linear functions exact
!-------------------------------------------------------------------     

MODULE matrixcorr
 IMPLICIT NONE
 REAL, DIMENSION(:,:,:), ALLOCATABLE :: gradmatrix
END MODULE matrixcorr

!-------------------------------------------------------------------
!  quantities used for link-list neighbour finding
!-------------------------------------------------------------------

MODULE linklist
 USE dimen_mhd
 IMPLICIT NONE
 INTEGER, DIMENSION(:), ALLOCATABLE :: ll,ifirstincell,iamincell,numneigh
 INTEGER, DIMENSION(ndim) :: ncellsx
 INTEGER :: ncells,ncellsloop
 REAL :: dxcell
END MODULE linklist

!-------------------------------------------------------------------
!  logical unit numbers for input and output files
!-------------------------------------------------------------------

MODULE loguns
 IMPLICIT NONE
 INTEGER :: iprint,ievfile,idatfile,iread,ireadf
 INTEGER :: ifile
 CHARACTER(LEN=120) :: rootname       ! name of the run
END MODULE loguns   

!-------------------------------------------------------------------
!  program options
!-------------------------------------------------------------------

MODULE options
 use dimen_mhd
 IMPLICIT NONE
 INTEGER :: iener,icty,iav,ikernav
 INTEGER :: iprterm,idumpghost,ihvar
 INTEGER :: imhd,imagforce,idivBzero   !  (mhd options)
 INTEGER :: iexternal_force,ixsph,isplitpart
 INTEGER :: igravity,ikernel,ikernelalt,iresist
 INTEGER :: maxdensits
 INTEGER, DIMENSION(ndim) :: ibound
 INTEGER, DIMENSION(3) :: iavlim
 REAL :: damp,dampz,dampr,psidecayfact,tolh,hsoft,etamhd,rhocrit
 CHARACTER(LEN=12) :: geom
 LOGICAL :: usenumdens
END MODULE

!-------------------------------------------------------------------
!  basic particle properties
!-------------------------------------------------------------------

MODULE part
 USE dimen_mhd
 IMPLICIT NONE
 INTEGER :: npart,nbpts,ntotal
 INTEGER, DIMENSION(:), ALLOCATABLE :: itype
 REAL, DIMENSION(:), ALLOCATABLE :: pmass,sqrtg
 REAL, DIMENSION(:,:), ALLOCATABLE :: x   
 REAL, DIMENSION(:), ALLOCATABLE :: dens,rho,pr,uu,en,hh,psi,spsound
 REAL, DIMENSION(:,:), ALLOCATABLE :: vel,pmom,sourceterms,alpha
 REAL, DIMENSION(:,:), ALLOCATABLE :: Bfield, Bevol
 REAL, DIMENSION(ndimB) :: Bconst
END MODULE part 

!-------------------------------------------------------------------
!  particle properties at the beginning of the time step
!-------------------------------------------------------------------

MODULE part_in
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: rhoin,prin,hhin,enin,psiin
 REAL, DIMENSION(:,:), ALLOCATABLE :: xin,velin,pmomin,alphain
 REAL, DIMENSION(:,:), ALLOCATABLE :: Bevolin
END MODULE

!-------------------------------------------------------------------
!  rates of change of particle properties (force, derivatives)
!-------------------------------------------------------------------

MODULE rates
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: drhodt,dudt,dendt,dhdt,dpsidt,poten
 REAL, DIMENSION(:,:), ALLOCATABLE :: force,dBevoldt,daldt,gradpsi
 REAL :: potengrav
END MODULE rates

!-------------------------------------------------------------------
!  rates of change at the beginning of the time step
!-------------------------------------------------------------------

!MODULE rates_in      ! only needed if using leapfrog
! IMPLICIT NONE
! REAL, DIMENSION(:), ALLOCATABLE :: drhodtin,dudtin,dendtin,daldtin
! REAL, DIMENSION(:,:), ALLOCATABLE :: forcein,dBfielddtin
!END MODULE rates_in

!-------------------------------------------------------------------
!  initial particle separation, initial smoothing length (in units of psep)
!-------------------------------------------------------------------

MODULE setup_params
 USE dimen_mhd
 USE options, only:geom
 IMPLICIT NONE
 CHARACTER(LEN=12) :: geomsetup
 REAL, PARAMETER :: pi = 3.1415926536
!
!--parameters for 2D-MRI simulations
!
 REAL, PARAMETER :: Rcentre = 100.
 REAL, PARAMETER :: Omega2 = 1./Rcentre**3
 REAL, PARAMETER :: Omega0 = Rcentre**(-1.5)
 REAL, PARAMETER :: domegadr = 1.5
 REAL :: psep,hfact, R_grav,xlayer,Alayercs,dwidthlayer
 real, parameter :: omegafixed = 1.0
 
END MODULE setup_params

!-------------------------------------------------------------------
!  time stepping related quantities
!-------------------------------------------------------------------

MODULE timestep
 IMPLICIT NONE
!
!  (max tsteps, number of steps before output, current number of steps,
!   number of steps before using a direct summation)   
!
 INTEGER :: nmax,nout,nsteps,ndirect,nsubsteps_divB,iseedMC
!
! (max time, time before output, dt)
!
 REAL :: tmax,tout,dt,dt0,time,vsig2max,dtscale,w0,Omega0
!
! (time step criterion from forces, courant condition)
!
 REAL :: dtforce, dtcourant, C_force, C_cour, C_rho, dtrho
 REAL :: dtav
END MODULE timestep

!-------------------------------------------------------------------
!   version number
!-------------------------------------------------------------------

MODULE versn
 IMPLICIT NONE
 CHARACTER(LEN=30) :: version
END MODULE versn

!-------------------------------------------------------------------
! XSPH factor
!-------------------------------------------------------------------

MODULE xsph
 IMPLICIT NONE
 REAL, DIMENSION(:,:), ALLOCATABLE :: xsphterm
 REAL :: xsphfac
END MODULE xsph
