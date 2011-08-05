!----------------------------------------------------------------
!     Set up a uniform density grid of particles in 1D
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
 INTEGER :: imax
 REAL :: sigma,massp,totmass,pri
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
 IF (ndim .GT. 1) STOP 'unifdis not implemented for ndim > 1'
!
!--initially set up a uniform density grid
! 	    
 ibound = 1	! boundary type
 nbpts = 6	! use ghosts not fixed
 xmin(1) = -0.5	! set position of boundaries
 xmax(1) = 0.5
 imax = INT((xmax(1)-xmin(1))/psep)
 
 totmass = 1.0	!4./3.
 massp = totmass/imax	! average particle mass
 sigma = 1./SQRT(2.)
 iexternal_force = 0
!
!--allocate memory here
!
 CALL alloc(imax)
 
 DO i=1,imax
    x(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    vel(:,i) = 0.
    IF (x(1,i).LT.0.) THEN
       dens(i) = 10.0
       pmass(i) = 10.0*massp
    ELSE
       dens(i) = 1.0
       pmass(i) = massp
    ENDIF
    pri = 1.0
    uu(i) = pri/((gamma-1.)*dens(i))
    IF (imhd.NE.0) THEN 
       Bfield(1,i) = 0.
       Bfield(2,i) = sigma*dens(i)
       Bfield(3,i) = 0.
    ENDIF 
!   print*,i,x(i),dens(i),uu(i),pmass(i),x(i)-x(i-1)
 ENDDO
  
 npart = imax
 ntotal = npart
 
 RETURN
END
