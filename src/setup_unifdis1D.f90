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
 REAL :: sigma,massp,totmass
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
 IF (ndim .GT. 1) STOP 'unifdis not implemented for ndim > 1'
!
!--initially set up a uniform density grid
! 	    
 ibound = 0	! periodic boundaries
 nbpts = 0		! use ghosts not fixed
 xmin(1) = -1.0	! set position of boundaries
 xmax(1) = 1.0
 imax = INT((xmax(1)-xmin(1))/psep)
 
 totmass = 4./3.
 massp = totmass/imax	! average particle mass
 sigma = 1./SQRT(2.)
 itoystar = 1
!
!--allocate memory here
!
 CALL alloc(imax,1)
 
 DO i=1,imax
    xin(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    velin(:,i) = 0.
    rhoin(i) = 1.0
    pmass(i) = massp
    uuin(i) = 0.3
    hhin(i) = hfact*(massp/rhoin(i))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
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
