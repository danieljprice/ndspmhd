!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the Bx peak advection problem in Dedner et al JCP 175, 645  !!
!!  Bx = r(x^2 + y^2)/sqrt(4pi) (ie div B .ne. 0)                         !!
!!  Basically to see how an initially non-zero div B propagates           !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE eos
 USE options
 USE part
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,npartx,nparty,ipart
 REAL :: denszero,przero
 REAL :: pri,rbump,rr
 REAL :: totmass,gam1,massp,const
 REAL, DIMENSION(ndim) :: xorigin, dx
 REAL, DIMENSION(ndimV) :: Bzero
!
!--check number of dimensions is right
!
 IF (ndim.NE.2) STOP ' ndim must be = 2 for MHD rotor problem'
 IF (ndimV.NE.3) STOP ' need ndimV=3 for this problem'
!
!--set boundaries
!            	    
 ibound = 3     ! reflective ghosts (boundaries not important in this problem)
 nbpts = 0      ! no fixed particles
 xmin(:) = -0.5 ! unit square
 xmax(:) = 1.5
 const = SQRT(4.*pi) 
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0 ! co-ordinates of the centre of the initial blast
 rbump = 0.125		! radius of the initial bump
 Bzero(:) = 0.
 IF (imhd.NE.0) Bzero(3) = 1.0/const	! uniform field in Bz direction
 przero = 6.0		! initial pressure
 denszero = 1.0		! ambient density
 
 gam1 = gamma - 1.

 WRITE(iprint,*) 'Two dimensional div B advection problem '
 WRITE(iprint,10) denszero,rbump,Bzero(3),przero
10 FORMAT(/,' density  = ',f10.3,', size of bump = ',f6.3,/, &
            ' Initial Bz   = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,psep,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass in ambient medium
!
 totmass = denszero*PRODUCT(xmax(:)-xmin(:))
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    dx(:) = x(:,ipart)-xorigin(:) 
    rr = DOT_PRODUCT(dx,dx)
    Bfield(:,ipart) = Bzero(:)
    IF (rr.LE.rbump) THEN
       Bfield(1,ipart) = (4096.*rr**4 - 128.*rr**2 + 1.)/const
    ENDIF  
    pmass(ipart) = massp
    dens(ipart) = denszero
    vel(:,ipart) = 0.
    vel(1:ndim,ipart) = 1.0
    pri = przero 
    uu(ipart) = pri/(gam1*denszero)
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
