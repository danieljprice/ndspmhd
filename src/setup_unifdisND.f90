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
 USE options
 USE part
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER :: i
 REAL :: massp,volume,totmass
 REAL :: denszero
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 3	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5
!
!--set up the uniform density grid
!
 CALL set_uniform_cartesian(4,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = PRODUCT(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    Bfield(:,i) = 0.
 ENDDO 
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
  
 RETURN
END
