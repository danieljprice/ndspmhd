!!-----------------------------------------------------------------
!!  Sets ghost particles if needed (ND)
!!  This version just loops over all the particles, and copies
!!  those which lie near the boundaries.
!!
!!  (previous version tried to use the link list to find particles
!!   near the boundary then created new link list cells)
!!
!! works in 2D, 3D need to make dxprev etc store previous 2 boundaries
!! for corners
!!
!!-----------------------------------------------------------------

SUBROUTINE set_ghost_particles
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE derivB
! USE linklist
 USE options
 USE part
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,k,kstart
 INTEGER :: ipart,jpart,jprev,idimen,idimenprev,ighost
 REAL :: tol,dxbound
 REAL :: dx,xbound,xperbound
 INTEGER, PARAMETER :: maxghost = 2*ndim	! max ghosts for any particle
 INTEGER, DIMENSION(ndim) :: nghost
 REAL, DIMENSION(maxghost) :: dxprev,xboundprev,xperboundprev
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine set_ghost_particles'
!
!--our link list grid used for neighbour finding bins particles within 2h
!  of the boundary in the first and last cells in the x direction
!      
 ntotal = npart	! reset to no ghost particles
!
!--use maximum value of h
! 
 hhmax = MAXVAL(hh(1:npart))
 dxbound = radkern*hhmax	! this is the maximum distance away from boundary
! PRINT*,' hhmax, hhmin = ',hhmax,MAXLOC(hh(1:npart)),MINVAL(hh(1:npart)),MINLOC(hh(1:npart))
 IF (ANY ((xmax(:)-xmin(:)) < dxbound)) THEN
    WRITE(iprint,*) ' hhmax too large = ',hhmax,MAXLOC(hh(1:npart))
    CALL quit
 ENDIF
!
!--loop over all the particles
!
 over_part: DO jpart=1,npart
    
    ighost = 0			! this particle initially has no ghosts
!
!--if using reflecting ghosts, reflect if dx < 2*hi, as particle sees own ghost
!  for periodic must use hmax as particle sees different particles as ghosts
!  (for efficiency should use max h of all particles near boundary)
!
    IF (ibound.EQ.2) dxbound = radkern*hh(jpart)
    
    over_dimen: DO idimen = 1, ndim		! over spatial dimensions
      
       nghost(idimen) = 0
      
       over_minmax: DO i = 1, 2			! over xmin, xmax  

          IF (i.EQ.1) THEN
	     xbound = xmin(idimen)
             xperbound = xmax(idimen)	! boundary where periodic ghosts go       
             dx = x(idimen,jpart) - xmin(idimen)
          ELSE
             xbound = xmax(idimen)
             xperbound = xmin(idimen)
	     dx = xmax(idimen) - x(idimen,jpart)
          ENDIF                 
!
!--copy if < 2h from boundary, don't copy if on or over the boundary (dx <= 0)
!  
          IF ((dx.LT.dxbound).AND.(dx.GT.0)) THEN 
!
!--count number of reflections in this dimension
!
	     nghost(idimen) = nghost(idimen) + 1
!
!--dx is now the amount to shift the particle from xbound/xperbound
!	     
	     dx = x(idimen,jpart) - xbound
	     
	     ighost = ighost + 1	! number of ghosts this particle has
!
!--save the boundary in case the particle is near an edge
!     	     
	     dxprev(ighost) = dx
             xboundprev(ighost) = xbound
	     xperboundprev(ighost) = xperbound
!	     idimenprev(ighost) = idimen
	     
	     ipart = ntotal + 1 	! label of new boundary particle
             IF (ipart.GT.SIZE(rho)) THEN
	        WRITE(iprint,*) 'ghost: ntotal > array size, re-allocating... '
	        CALL alloc(SIZE(rho),2)
	     ENDIF 
	     ntotal = ipart		! add one to total number of particles

             x(:,ipart) = x(:,jpart)		! copy positions /vels
	     vel(:,ipart) = vel(:,jpart)	! then overwrite
	     
             IF (ibound.EQ.2) THEN 		! reflecting
	        x(idimen,ipart) = xbound - dx
		vel(:,ipart) = -vel(:,jpart)	! should reflect all vels
	     ELSEIF (ibound.EQ.3) THEN		! periodic
	        x(idimen,ipart) = xperbound + dx
	     ENDIF
	     ireal(ipart) = jpart

!	     PRINT*,'copying  old particle x(',jpart,'),vel,rho =',	&
!	            x(:,jpart),vel(:,jpart),rho(jpart)
!	     PRINT*,'creating new particle x(',ipart,'),vel,rho =',	&
!	            x(:,ipart),vel(:,ipart),rho(ipart)
!
!--make additional ghost(s) if near an edge
!	     
             IF (idimen.GT.1) THEN
	      kstart = 1
	      DO idimenprev = 1,ndim-1
	       IF (nghost(idimenprev).GE.1) THEN	! if already ghosts in x = edge		

		DO k=kstart,nghost(idimenprev)
		
		  ipart = ntotal + 1
                  IF (ipart.GT.SIZE(rho)) THEN
	             WRITE(iprint,*) 'ghost: ntotal > array size, re-allocating... '
	             CALL alloc(SIZE(rho),2)
	          ENDIF 
                  ntotal = ipart
		
                  x(:,ipart) = x(:,jpart)		! copy positions /vels
	          vel(:,ipart) = vel(:,jpart)		! then overwrite
		  IF (ibound.EQ.2) THEN 		! reflecting
	             x(idimen,ipart) = xbound - dx
		     x(idimenprev,ipart) = xboundprev(k) - dxprev(k)
		     vel(:,ipart) = -vel(:,jpart)	! should reflect all vels
	          ELSEIF (ibound.EQ.3) THEN		! periodic
	             x(idimen,ipart) = xperbound + dx
		     x(idimenprev,ipart) = xperboundprev(k) + dxprev(k)
	          ENDIF
		  ireal(ipart) = jpart

!	          PRINT*,'edge: copying  old particle x(',jpart,'),vel,rho =',	&
!	              x(:,jpart),vel(:,jpart),rho(jpart),idimen
!	          PRINT*,'edge: creating new particle x(',ipart,'),vel,rho =',	&
!	              x(:,ipart),vel(:,ipart),rho(ipart)

                ENDDO
		
		kstart = nghost(idimenprev)
		
	       ENDIF	
	      ENDDO
	     ENDIF

!             xboundprev = xbound	! save the boundary in case it is
!	     xperboundprev = xperbound  !  near an edge
!	     dxprev = dx
!	     idimenprev = idimen
		    
          ENDIF
       ENDDO over_minmax
    ENDDO over_dimen
 ENDDO over_part
!
!--copy particle quantities to the ghost particles
!
 DO i=npart+1,ntotal
    j = ireal(i)
    pmass(i) = pmass(j)
    rho(i) = rho(j)
    uu(i) = uu(j)
    en(i) = en(j)
    hh(i) = hh(j)
    alpha(i) = alpha(j)
    Bfield(:,i) = Bfield(:,j)
    divB(i) = divB(j)
 ENDDO
!
!--set unused elements of the array to zero (can cause errors in eos)
! 
 DO i = ntotal+1,SIZE(rho)
    rho(i) = 0.
    uu(i) = 0.
 ENDDO

 RETURN
END SUBROUTINE set_ghost_particles
      
