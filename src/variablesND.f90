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
 REAL :: hhmax    
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
!  equation of state related quantities
!-------------------------------------------------------------------

MODULE eos
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: spsound
 REAL :: gamma
END MODULE

!-------------------------------------------------------------------
!  Lorentz force
!-------------------------------------------------------------------

MODULE fmagarray
 IMPLICIT NONE
 REAL, DIMENSION(:,:), ALLOCATABLE :: fmag
END MODULE

!-------------------------------------------------------------------
! Gravitational potential and force
!-------------------------------------------------------------------
MODULE gravity
 IMPLICIT NONE
 REAL, DIMENSION(:,:), ALLOCATABLE :: fgrav
END MODULE gravity

!-------------------------------------------------------------------
!  correction terms when using a spatially variable smoothing length
!-------------------------------------------------------------------     

MODULE hterms
 IMPLICIT NONE
 INTEGER :: itsdensity
 REAL, DIMENSION(:), ALLOCATABLE :: gradh, gradhaniso
END MODULE

!-------------------------------------------------------------------
!  kernel tables
!-------------------------------------------------------------------     

MODULE kernel
 IMPLICIT NONE
 INTEGER, PARAMETER :: ikern=4000    ! dimensions of kernel table
 REAL, DIMENSION(0:ikern) :: wij,grwij,wijaniso,grwijaniso
 REAL :: dq2table,ddq2table,radkern2,radkern
END MODULE

!
!--these tables used for plotting only (not in rates etc)
!

MODULE kernelextra
 USE kernel
 IMPLICIT NONE
 REAL, DIMENSION(0:ikern) :: grgrwij,grgrwijaniso
 CHARACTER(LEN=100) :: kernelname
END MODULE kernelextra

!-------------------------------------------------------------------
!  quantities used for link-list neighbour finding
!-------------------------------------------------------------------

MODULE linklist
 USE dimen_mhd
 IMPLICIT NONE
 INTEGER, DIMENSION(:), ALLOCATABLE :: ll,ifirstincell,iamincell
 INTEGER, DIMENSION(ndim) :: ncellsx
 INTEGER :: ncells,ncellsloop,nlistdim
 REAL :: dxcell
END MODULE linklist

!-------------------------------------------------------------------
!  logical unit numbers for input and output files
!-------------------------------------------------------------------

MODULE loguns
 IMPLICIT NONE
 INTEGER :: iprint,ievfile,idatfile,iread,ireadf     
 CHARACTER(LEN=20) :: rootname       ! name of the run
END MODULE loguns   

!-------------------------------------------------------------------
!  program options
!-------------------------------------------------------------------

MODULE options
 use dimen_mhd
 IMPLICIT NONE
 INTEGER :: iener,icty,iav,ikernav
 INTEGER :: iprterm,idumpghost,ihvar
 INTEGER :: imhd,imagforce,idivBzero      !  (mhd options)
 INTEGER :: iexternal_force,ixsph,ianticlump
 INTEGER :: igravity,ikernel
 INTEGER :: igeom,maxdensits,iunity
 INTEGER, DIMENSION(ndim) :: ibound
 INTEGER, DIMENSION(3) :: iavlim
 REAL :: damp,psidecayfact
END MODULE

!-------------------------------------------------------------------
!  basic particle properties
!-------------------------------------------------------------------

MODULE part
 IMPLICIT NONE
 INTEGER :: npart,nbpts,ntotal
 INTEGER, DIMENSION(:), ALLOCATABLE :: itype
 REAL, DIMENSION(:), ALLOCATABLE :: pmass,sqrtg
 REAL, DIMENSION(:,:), ALLOCATABLE :: x   
 REAL, DIMENSION(:), ALLOCATABLE :: dens,rho,pr,uu,en,hh,psi
 REAL, DIMENSION(:,:), ALLOCATABLE :: vel,pmom,sourceterms,alpha
 REAL, DIMENSION(:,:), ALLOCATABLE :: Bfield, Bevol
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
!  constant used for a polytropic/isothermal equation of state
!-------------------------------------------------------------------

MODULE polyconst
 IMPLICIT NONE
 REAL :: polyk
END MODULE

!-------------------------------------------------------------------
!  rates of change of particle properties (force, derivatives)
!-------------------------------------------------------------------

MODULE rates
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: drhodt,dudt,dendt,dhdt,dpsidt
 REAL, DIMENSION(:,:), ALLOCATABLE :: force,dBevoldt,daldt
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
 IMPLICIT NONE
 REAL, PARAMETER :: pi = 3.1415926536
 REAL :: psep,hfact,hpower   ! hpower set in initialise
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
 INTEGER :: nmax,nout,nsteps,ndirect
!
! (max time, time before output, dt)
!
 REAL :: tmax,tout,dt,dt0,time     
!
! (time step criterion from forces, courant condition)
!
 REAL :: dtforce, dtcourant, C_force, C_cour
END MODULE timestep

!-------------------------------------------------------------------
!   unity function and its derivative
!-------------------------------------------------------------------

MODULE unityfunc
 IMPLICIT NONE
 REAL, DIMENSION(:), ALLOCATABLE :: unity
 REAL, DIMENSION(:,:), ALLOCATABLE :: gradunity
END MODULE unityfunc


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

!-------------------------------------------------------------------
! anticlumping term settings
!-------------------------------------------------------------------

MODULE anticlumping
 IMPLICIT NONE
 INTEGER :: neps
 REAL :: eps
END MODULE anticlumping
