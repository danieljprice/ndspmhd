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
 USE part_in
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
 itoystar = 0
!
!--allocate memory here
!
 CALL alloc(imax)
 
 DO i=1,imax
    xin(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    velin(:,i) = 0.
    IF (xin(1,i).LT.0.) THEN
       rhoin(i) = 10.0
       pmass(i) = 10.0*massp
    ELSE
       rhoin(i) = 1.0
       pmass(i) = massp
    ENDIF
    pri = 1.0
    uuin(i) = pri/((gamma-1.)*rhoin(i))
    hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	 ! ie constant everywhere
    IF (imhd.NE.0) THEN 
       Bin(1,i) = 0.
       Bin(2,i) = sigma*rhoin(i)
       Bin(3,i) = 0.
    ENDIF 
!   print*,i,xin(i),rhoin(i),uuin(i),pmass(i),xin(i)-xin(i-1)
 ENDDO
  
 npart = imax
 ntotal = npart
 
 RETURN
END
