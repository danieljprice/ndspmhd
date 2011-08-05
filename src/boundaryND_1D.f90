!!-------------------------------------------------------------------------
!! update of particles at end of step for periodic (particles crossing
!! domain) and inflow/outflow (adds / subtracts particles from the domain)
!!
!! for ND case only periodic boundaries implemented
!!
!! Changes log:
!! 17/10/03 - bug fix in inflow boundaries - forces/drhodt etc copied
!!-------------------------------------------------------------------------
	 
SUBROUTINE boundary
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE fmagarray
 USE options
 USE part
 USE part_in
 USE rates
 USE xsph
! USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,nnew,nsub,nsubtemp,jdim
 REAL :: psepleft,psepright
 LOGICAL :: debugging
!
!--allow for tracing flow
!      
 IF (trace) THEN
    WRITE(iprint,*) ' Entering subroutine boundary'
    debugging = .true.
 ENDIF   
 
!---------------------------------------------------------------------------
!  inflow/outflow when fixed particles are used.
!---------------------------------------------------------------------------     
 IF (ibound.EQ.1) THEN

    IF (ndim.GT.1) THEN    
!       WRITE(iprint,*) 'warning: inflow/outflow not implemented in ND'
    ELSE
!
!--this is just pulled straight from the 1D code - needs to be
!  modified if applied in 2 or 3D since it assumes that particles
!  are ordered from left to right in the 1D domain

!
!--determine particle separation at left and right boundaries
!
    psepleft = x(1,2) - x(1,1)
    psepright = x(1,npart) - x(1,npart-1)
    
    IF (psepright.LT.0.) WRITE(iprint,*) 'Boundary: warning psepright < 0', &
                               psepright,x(1,npart-1),x(1,npart),npart
    IF (psepleft.LT.0.) WRITE(iprint,*) 'Boundary: warning psepleft < 0',psepleft,x(1,1),x(1,2)
!
!--inflow to left boundary
!
    IF (x(1,1).GT.(xmin(1)+psepleft)) THEN
       IF (debugging) WRITE(iprint,*) 'inflow to left boundary npart = ',npart+1
       IF ((x(1,1)-xmin(1)).GT.2*psepleft) THEN
          WRITE(iprint,*)' Need more than one particle',psepleft,(x(1,1)-xmin(1))/psepleft
	  nnew = INT((x(1,1)-xmin(1))/psepleft)
	  WRITE(iprint,*)' Number of particles needed = ',nnew
       ENDIF	  
       DO i=npart,1,-1		! relabel particles
!          print*,' particle ',i+1,' = ',i
	  x(1,i+1) = x(1,i)		! (must copy all particle properties)
	  vel(:,i+1) = vel(:,i)
	  rho(i+1) = rho(i)
	  hh(i+1) = hh(i)
	  uu(i+1) = uu(i)
	  en(i+1) = en(i)
	  Bfield(:,i+1) = Bfield(:,i)
	  alpha(i+1) = alpha(i)
	  pmass(i+1) = pmass(i)
        
	  force(:,i+1) = force(:,i)
          drhodt(i+1) = drhodt(i)
          dudt(i+1) = dudt(i)
          dendt(i+1) = dendt(i)
          dBfielddt(:,i+1) = dBfielddt(:,i)
	  dhdt(i+1) = dhdt(i)
	  daldt(i+1) = daldt(i)
          xsphterm(:,i+1) = xsphterm(:,i)	! after here not crucial
          fmag(:,i+1) = fmag(:,i)
          divB(i+1) = divB(i)
          curlB(:,i+1) = curlB(:,i)
       ENDDO       
!
!--now make new particle number 1
!  (particle properties are automatically set in the next iteration of step)
!
       npart = npart + 1	! add new particle
       x(1,1) = x(1,1) - psepleft
       vel(:,1) = vel(:,2)	! copy quantities from particle 2
       rho(1) = rho(2)
       hh(1) = hh(2)
       uu(1) = uu(2)
       en(1) = en(2)
       Bfield(:,1) = Bfield(:,2)
       alpha(1) = alpha(2)
       pmass(1) = pmass(2)
       
!       PRINT*,' new particle x(',1,') = ',x(1)
!
!--outflow from left boundary
!       
    ELSEIF (x(1,1).LT.xmin(1)) THEN
       nsub = INT((xmin(1)-x(1,1))/psepleft) !+ 1
       IF (debugging) 	&
          WRITE(iprint,*) 'outflow from left boundary npart = ',npart-nsub,nsub
       npart = npart - nsub	! subtract particle(s)
       DO i=1,npart		! relabel particles
          x(1,i) = x(1,i+nsub)		! (must copy all particle properties)
	  vel(:,i) = vel(:,i+nsub)
	  rho(i) = rho(i+nsub)
	  hh(i) = hh(i+nsub)
	  uu(i) = uu(i+nsub)
	  en(i) = en(i+nsub)
	  Bfield(:,i) = Bfield(:,i+nsub)
	  alpha(i) = alpha(i+nsub)
	  pmass(i) = pmass(i+nsub)

	  force(:,i) = force(:,i+nsub)
          drhodt(i) = drhodt(i+nsub)
          dudt(i) = dudt(i+nsub)
          dendt(i) = dendt(i+nsub)	  
          dBfielddt(:,i) = dBfielddt(:,i+nsub)
	  dhdt(i) = dhdt(i+nsub)
	  daldt(i) = daldt(i+nsub)
          xsphterm(:,i) = xsphterm(:,i+nsub)	! after here not crucial
          fmag(:,i) = fmag(:,i+nsub)
          divB(i) = divB(i+nsub)
          curlB(:,i) = curlB(:,i+nsub)
       ENDDO
              
    ENDIF
!
!--right boundary
!    
    IF (x(1,npart).GT.xmax(1)) THEN		! outflow from right boundary
       nsub = INT((x(1,npart)-xmax(1))/psepright) !+ 1
       IF (debugging) WRITE(iprint,*) 'outflow from right boundary npart = ',npart-nsub,nsub       
!--if more than one particle, check other particles close to boundary
       IF (nsub.GE.2) THEN
          WRITE(iprint,*) 'need more than one particle ',nsub
          nsubtemp = 0
	  DO i=0,nsub
	     IF (x(1,npart-i).GT.xmax(1)) THEN
	        nsubtemp = nsubtemp + 1
	     ENDIF
	  ENDDO
	  nsub = nsubtemp
          WRITE(iprint,*) 'really need more than one particle ',nsub
       ENDIF	  
       npart = npart - nsub	! subtract particle(s)
    ELSEIF (x(1,npart).LT.(xmax(1)-psepright)) THEN	! inflow from right boundary   
       IF (debugging) WRITE(iprint,*) 'inflow to right boundary npart = ',npart+1
       IF ((x(1,npart)-xmax(1)).GT.2*psepright) THEN
          WRITE(iprint,*)' Need more than one particle',psepright,(x(1,npart)-xmax(1))/psepright
	  nnew = INT((x(1,npart)-xmax(1))/psepright)
	  WRITE(iprint,*)' Number of particles needed = ',nnew
       ENDIF	  

       npart = npart + 1	! add only one new particle
!
!--create new particle number npart
!  (force etc calculated next step - doesn't matter so long as its a fixed
!   particle to start with)
!
       x(1,npart) = x(1,npart-1) + psepright
       pmass(npart) = pmass(npart-1)
       rho(npart) = rho(npart-1)
       vel(:,npart) = vel(:,npart-1)
       hh(npart) = hh(npart-1)
       uu(npart) = uu(npart-1)
       en(npart) = en(npart-1)
       Bfield(:,npart) = Bfield(:,npart-1)
       alpha(npart) = alpha(npart-1)         
    ENDIF
 
    ntotal = npart
!    WRITE(iprint,*) ' New number of particles = ',npart
!
!--make sure all the initial quantities are the same
!
    DO i=1,npart
       xin(1,i) = x(1,i)
       velin(:,i) = vel(:,i)
       rhoin(i) = rho(i)
       hhin(i) = hh(i)
       uuin(i) = uu(i)
       enin(i) = en(i)
       Bfieldin(:,i) = Bfield(:,i)
       alphain(i) = alpha(i)
!
!--check to see if particles are ordered left to right
!       
!       IF (i.LT.npart.AND.x(1,i).GT.x(1,i+1)) THEN
!          WRITE(iprint,*) 'Warning: particles crossed ',i,x(1,i),i+1,x(1,i+1)       
!       ENDIF
    ENDDO
        
    ENDIF

!----------------------------------------------------------------------
! periodic boundary conditions - allow particles to cross the domain
!----------------------------------------------------------------------

 ELSEIF (ibound.EQ.3) THEN
    
    DO i=1,npart
! -- this is very f77 but I can''t get it right using where statements 
       DO jdim=1,ndim
          IF (x(jdim,i).GT.xmax(jdim)) THEN
!	     print*,'ss xold,xmax,xnew = ',jdim,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
	     x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
	     xin(jdim,i) = x(jdim,i)
	  ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!	     print*,'ss xold,xmin,xnew = ',jdim,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)	     
!	     read*
             x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
	     xin(jdim,i) = x(jdim,i)
          ENDIF	  
       ENDDO
    
    ENDDO   

 ENDIF

 RETURN 
END SUBROUTINE boundary
