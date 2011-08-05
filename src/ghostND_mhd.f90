!!-----------------------------------------------------------------
!!  Sets ghost particles if needed (ND)
!!  This version just loops over all the particles, and copies
!!  those which lie near the boundaries.
!!
!!  (previous version tried to use the link list to find particles
!!   near the boundary then created new link list cells)
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
 INTEGER, PARAMETER :: maxbound = 2	! maximum no. of boundaries in each dim.
 INTEGER, DIMENSION(ndim) :: nbound	! number of boundaries in each dim.
 INTEGER :: i,j,ighostprev,istart,ighostx
 INTEGER :: imaxmin,imaxminprev,imaxminprevprev
 INTEGER :: idimen,idimenprev,idimenprevprev
 INTEGER :: ipart,jpart,jprev,jtemp
 REAL, DIMENSION(ndim) :: dxbound,xpart
 REAL, DIMENSION(ndim,maxbound) :: xnew
 REAL :: dx,dxshift,xbound,xperbound
 LOGICAL, DIMENSION(ndim,maxbound) :: imakeghost
 LOGICAL :: ireflect
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine set_ghost_particles'
!
!--reset to zero ghost particles
!      
 ntotal = npart
 nbound(:) = 2 		! 2 boundaries in each dimension (ie. cartesian)
			! this could be changed elsewhere
			! must not be > maxbound
 jtemp = 1861
! PRINT*,'jtemp = ',jtemp
!
!--use maximum value of h - check 2h is not bigger than box size
! 
 hhmax = MAXVAL(hh(1:npart))
 dxbound(:) = radkern*hhmax	! this is the maximum distance away from boundary
! PRINT*,' hhmax, hhmin = ',hhmax,MAXLOC(hh(1:npart)),MINVAL(hh(1:npart)),MINLOC(hh(1:npart))
 IF (ANY ((xmax(:)-xmin(:)) < dxbound)) THEN
    WRITE(iprint,*) 'WARNING: ghost: hhmax too large = ',hhmax,MAXLOC(hh(1:npart))
 ENDIF
!
!--loop over all the particles (with link list could supply a list of 
!  particles within 2h of the boundary to search)
!
 over_part: DO jpart=1,npart
!
!--if using reflecting ghosts, reflect if dx < 2*hi, as particle sees own ghost
!  for periodic must use hmax as particle sees different particles as ghosts
!  (for efficiency in this case should use max h of all particles near boundary)
!
    WHERE (ibound.EQ.2) dxbound(:) = radkern*hh(jpart)
!    IF (jpart.EQ.jtemp) PRINT*,jpart,' x = ',x(:,jpart)
    
    over_dimen: DO idimen = 1, ndim		! over spatial dimensions

       ireflect = .false.
      
       over_maxmin: DO imaxmin = 1, nbound(idimen)	! over xmax, xmin  
!
!--set ghost position initially equal to particle position
!
          xpart(:) = x(:,jpart)
!
!--compute distance from boundary (dx)
!
	  IF (imaxmin.EQ.1) THEN	! max
             xbound = xmax(idimen)	! xbound is current boundary
             xperbound = xmin(idimen)	! boundary where periodic ghosts go
	     dx = xmax(idimen) - x(idimen,jpart)
          ELSE				! min
	     xbound = xmin(idimen)
             xperbound = xmax(idimen)
             dx = x(idimen,jpart) - xmin(idimen)
          ENDIF
!
!--copy if dx < 2h, don't copy if on or over the boundary (dx <= 0)
!  
          imakeghost(idimen,imaxmin) = ((dx.LT.dxbound(idimen)).AND.(dx.GT.0))
          
	  IF (imakeghost(idimen,imaxmin)) THEN
!
!--dxshift is now the amount to shift the particle from xbound/xperbound
!  (zero shift in other dimensions)
!	     
	     dxshift = x(idimen,jpart) - xbound
!
!--xnew is the shifted position of the ghost particle
!  save this for each boundary to use for edges/corners
!	     
	     IF (ibound(idimen).EQ.3) THEN	! periodic
	        xnew(idimen,imaxmin) = xperbound + dxshift
	     ELSEIF (ibound(idimen).EQ.2) THEN	! reflective
	        ireflect = .true.
	        xnew(idimen,imaxmin) = xbound - dxshift
	     ENDIF
!
!--set ghost position in current dimension equal to shifted position
!	     
	     xpart(idimen) = xnew(idimen,imaxmin)
!	     IF (jpart.EQ.jtemp) THEN
!	        PRINT*,jpart,' dim=',idimen,'bnd=',imaxmin,' ghost = ',xpart
!	     ENDIF
!
!--make the ghost particle at this position
!
	     CALL makeghost(jpart,xpart,ireflect)
!
!--make additional ghost(s) if near an edge/corner	     
!
! loop over previous dimensions
! (so if currently doing y, check y-x; if doing z, check z-y, z-x)
!
             IF (idimen.GT.1) THEN	! only in 2D or 3D
	      DO idimenprev = 1,idimen-1	! over previous dimension(s)
	        DO imaxminprev=1,nbound(idimenprev)	! over boundaries
!
!--make an edge ghost if there are ghosts set in this previous dimension
! (could be up to 12 edge reflections of this particle in 3D)
!	       
		   IF (imakeghost(idimenprev,imaxminprev)) THEN
!
!--set ghost particle position in previous dimension to shifted position
!		      
		      xpart(idimenprev) = xnew(idimenprev,imaxminprev)
!		      IF (jpart.EQ.jtemp) THEN
!		         PRINT*,'edge', idimen,idimenprev,'bnd=',imaxmin,imaxminprev,xpart
!		         PRINT*,'i_ghost = ',ntotal+1
!		      ENDIF
!
!--make the ghost particle at this position
!
		      IF (ibound(idimenprev).EQ.2) ireflect = .true.
		      CALL makeghost(jpart,xpart,ireflect)
!
!--make a corner ghost if ghosts set in both previous dimensions
!  (could be up to 8 corner reflections of this particle in 3D
!               
	              IF (idimenprev.GE.2) THEN
		        idimenprevprev = idimenprev-1	! in 3D
		        DO imaxminprevprev = 1,nbound(idimenprevprev)
			   IF (imakeghost(idimenprevprev,imaxminprevprev)) THEN
! 
!--set ghost particle position in previous dimension to shifted position
!		          
		              xpart(idimenprevprev) = xnew(idimenprevprev,imaxminprevprev)
!
!--make the ghost particle at this position
!
!			      IF (jpart.EQ.jtemp) THEN
!			         PRINT*,'corner',idimen,idimenprev,idimenprevprev				         
!			      ENDIF	 
			      IF (ibound(idimenprevprev).EQ.2) ireflect = .true.
		              CALL makeghost(jpart,xpart,ireflect)
!--reset boundary
			      xpart(idimenprevprev) = x(idimenprevprev,jpart)
			   ENDIF	! made ghost in previous dimension
			ENDDO	! over boundaries
		      ENDIF	! if > 2D
!
!--reset position in idimenprev
!			      
		      xpart(idimenprev) = x(idimenprev,jpart)			   

		      
		   ENDIF   ! made ghost in previous dimension
	        ENDDO	! over boundaries

	      ENDDO	! over previous dimensions	     
	     ENDIF	! if > 1D
	      
		    
          ENDIF	! make ghost for this boundary
       ENDDO over_maxmin
    ENDDO over_dimen
   
 ENDDO over_part
!
!--copy particle quantities to the ghost particles (conservative only)
!  (nb: the conservative quantities have not been set in the first call to ghosts)
!
 DO i=npart+1,ntotal
    j = ireal(i)
    pmass(i) = pmass(j)
    rho(i) = rho(j)
    hh(i) = hh(j)
    alpha(i) = alpha(j)
    psi(i) = psi(j)

    en(i) = en(j)
    Bcons(:,i) = Bcons(:,j)
    divB(i) = divB(j)
!    IF (j.EQ.jtemp) PRINT*,' ghost ',i
 ENDDO
!
!--set unused elements of the array to zero (can cause errors in eos)
! 
 IF (SIZE(rho).GT.ntotal) THEN
    DO i = ntotal+1,SIZE(rho)
       rho(i) = 0.
       uu(i) = 0.
    ENDDO
 ENDIF

 RETURN
END SUBROUTINE set_ghost_particles

!
! subroutine to make a ghost particle
!
SUBROUTINE makeghost(jpart,xghost,ireflect)
 USE dimen_mhd
 USE bound
! USE derivB
 USE loguns
 USE part
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: jpart	! index of particle to be ghosted
 INTEGER :: ipart
 REAL, INTENT(IN), DIMENSION(ndim) :: xghost ! position of ghost particle
 LOGICAL, INTENT(IN) :: ireflect
!
!--create new particle and reallocate memory if needed
! 
 ipart = ntotal + 1
 IF (ipart.GT.SIZE(rho)) THEN
    WRITE(iprint,*) 'ghost: ntotal > array size, re-allocating... '
    CALL alloc(SIZE(rho))
 ENDIF 
 ntotal = ipart
!
!--set position of new ghost particle
!
 x(:,ipart) = xghost(:)
!
!--copy velocities
!
 IF (ireflect) THEN 	 		! reflecting
    vel(:,ipart) = -vel(:,jpart)	! should reflect all vels
 ELSE					! periodic/fixed
    vel(:,ipart) = vel(:,jpart)
 ENDIF
!
!--ireal for ghosts refers to the real particle of which they are ghosts
!  
 ireal(ipart) = jpart  
!
!--copy particle properties
!
! pmass(ipart) = pmass(jpart)
! rho(ipart) = rho(jpart)
! uu(ipart) = uu(jpart)
! en(ipart) = en(jpart)
! hh(ipart) = hh(jpart)
! alpha(ipart) = alpha(jpart)
! psi(ipart) = psi(jpart)
! Bcons(:,ipart) = Bcons(:,jpart)
! divB(ipart) = divB(jpart)

! PRINT*,'copying  old particle x(',jpart,'),vel,rho =',	&
!        x(:,jpart),vel(:,jpart),rho(jpart)
! PRINT*,'creating new particle x(',ipart,'),vel,rho =',	&
!        x(:,ipart),vel(:,ipart),rho(ipart)
 
 RETURN
END SUBROUTINE makeghost

      
