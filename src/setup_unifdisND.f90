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
 USE part_in
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER :: i
 REAL :: massp,volume,totmass
 REAL :: rhozero
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
 IF (ndim .GT. 1) STOP 'unifdis not implemented for ndim > 1'
!
!--set boundaries
! 	    
 ibound = 2	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5
!
!--set up the uniform density grid
!
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass
!
 rhozero = 1.0
 volume = PRODUCT(xmax(:)-xmin(:))
 totmass = rhozero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO i=1,ntotal
    velin(:,i) = 0.
    rhoin(i) = rhozero
    pmass(i) = massp
    uuin(i) = 1.0	! isothermal
    hhin(i) = hfact*(massp/rhoin(i))**hpower	 ! ie constant everywhere
    Bin(:,i) = 0.
 ENDDO 
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
  
 RETURN
END
