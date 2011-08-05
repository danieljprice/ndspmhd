!!-----------------------------------------------------------------
!!  Sets ghost particles if needed (ND)
!!-----------------------------------------------------------------

SUBROUTINE set_ghost_particles
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE linklist
 USE options
 USE part
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j
 INTEGER :: icell,ibcell,ipart,jpart,jprev 
 REAL :: hhmax,tol
 REAL, DIMENSION(ndim) :: dx,xbound,xperbound
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
!--loop over cells which are close to boundaries
!      
 over_cells: DO j = 1, ncellsx(2)	! number of boundary cells
    
    min_max: DO i = 1, 2		! over xmin, xmax
    
       IF (MOD(i,2).EQ.0) THEN
          icell = (j-1)*ncellsx(1) + 1
          IF (ibound.EQ.2) THEN	! reflective
	     ibcell = icell - 1		! put ghosts in cell 0
          ELSEIF (ibound.EQ.3) THEN	! periodic
	     ibcell = ncells + 1		! put ghosts in cell ncells+1
          ENDIF
          xbound(1) = xmin(1)
          xperbound(1) = xmax(1)	! boundary where periodic ghosts go
       
       ELSE
          icell = j*ncellsx(1)
          IF (ibound.EQ.2) THEN	! reflective
             ibcell = ncells + 1	! put ghosts in cell ncells + 1
          ELSEIF (ibound.EQ.3) THEN	! periodic
             ibcell = 0		! put ghosts in cell 0
          ENDIF
          xbound(1) = xmax(1)
          xperbound(1) = xmin(1)
       ENDIF      

       jpart = ifirstincell(icell)
       ifirstincell(ibcell) = -1		! initialise head of chain
       jprev = jpart

       IF (jpart.EQ.-1) THEN
          WRITE(iprint,*) ' Warning: no particles near boundary cell ',icell
       ELSE
          DO WHILE(ll(jprev).NE.-1)
             dx(:) = x(:,jpart)-xbound(:)
             IF (abs(dx(1)).GT.0.) THEN	! don't copy if on the boundary
	     ipart = ntotal + 1	! label of new boundary particle
	     ntotal = ipart		! ntot = ntot + 1
	     IF (ntotal.GT.idim) THEN
	        WRITE(iprint,*) 'ghost: ntotal>idim, exit'
	        CALL quit
             ENDIF
             ll(ipart) = ifirstincell(ibcell)	! link list to old head of chain
	     ifirstincell(ibcell) = ipart 	! move head of chain to current particle
            
	     IF (ibound.EQ.2) THEN 		! reflecting
	        x(:,ipart) = xbound(:) - dx(:)
	        vel(:,ipart) = -vel(:,jpart)	! should reflect all vels
	     ELSEIF (ibound.EQ.3) THEN	! periodic
	        x(:,ipart) = xperbound(:) + dx(:)
	        vel(:,ipart) = vel(:,jpart)
	     ENDIF
	       
	     pmass(ipart) = pmass(jpart)
	     rho(ipart) = rho(jpart)
	     uu(ipart) = uu(jpart)
	     en(ipart) = en(jpart)
	     hh(ipart) = hh(jpart)
	     alpha(ipart) = alpha(jpart)
	     IF (imhd.NE.0) THEN
	        Bfield(:,ipart) = Bfield(:,jpart)
	        divB(ipart) = divB(jpart)
	     ENDIF
	     ireal(ipart) = jpart	
          ELSE
             WRITE(iprint,*) 'ghost: particle on boundary x = ',x(:,jpart),jpart
          ENDIF
          jprev = jpart
          IF (jprev.NE.-1) jpart = ll(jpart)
      ENDDO
    ENDIF
      
    ENDDO min_max	! over xmin, xmax
 ENDDO over_cells ! over number of boundary cells
            
 RETURN
END SUBROUTINE set_ghost_particles
      
