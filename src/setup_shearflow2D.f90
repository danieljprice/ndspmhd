!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

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
 INTEGER :: i
 REAL :: massp,volume,totmass
 REAL :: denszero,uuzero,cs0,polyk0
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
 WRITE(iprint,*) '2D Cartesian shear flow'
!
!--set boundaries
! 	    
 ibound = 3	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = 0.0  ! set position of boundaries
 xmax(:) = 1.0
!
!--set up the uniform density grid
!
 CALL set_uniform_cartesian(11,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = PRODUCT(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--sound speed
!
 cs0 = 0.05
 uuzero = cs0**2/(gamma*(gamma-1))
 polyk0 = cs0**2/(gamma*denszero**(gamma-1.))
 WRITE(iprint,*) ' cs0 = ',cs0, ' u = ',uuzero
 WRITE(iprint,*) ' polyk = ',polyk,' should be = ',polyk0
!
!--now assign particle properties
! 
 DO i=1,ntotal
    vel(:,i) = 0.
    vel(2,i) = SIN(2.*pi*x(1,i))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero ! isothermal
    Bfield(:,i) = 0.
 ENDDO 
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
  
 RETURN
END
