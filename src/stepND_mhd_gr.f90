!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!--------------------------------------------------------------------
	 
SUBROUTINE step
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE eos
 USE hterms
 USE options
 USE part
 USE part_in
 USE rates
 USE timestep
 USE setup_params
 USE xsph
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,jdim,ikernavprev,ierr
 REAL :: hdt
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine step'
!
!--set initial quantities
!
 hdt = 0.5*dt
 DO i=1,npart
    xin(:,i) = x(:,i)
    pmomin(:,i) = pmom(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
 ENDDO   
!
!--Mid-point Predictor step
!      
 DO i=1,npart
    IF (itype(i).EQ.1) THEN	! fixed particles
       pmom(:,i) = pmomin(:,i)
       rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)	     
       IF (iener.NE.0) en(i) = enin(i)
       hh(i) = hhin(i)	    
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(:,i) = alphain(:,i)
    ELSE
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       pmom(:,i) = pmomin(:,i) + hdt*force(:,i)
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
       IF (ihvar.EQ.1) THEN
!	   hh(i) = hfact(pmass(i)/rho(i))**hpower	! my version
	  hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower		! Joe's	   
       ELSEIF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (iavlim.NE.0) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)	   
    ENDIF
!
!--for periodic boundaries, allow particles to cross the domain
!  (this is only temporary as it is for the predicted quantity)
!    
    IF (ANY(ibound.EQ.3)) THEN
       DO jdim=1,ndim
          IF (ibound(jdim).EQ.3) THEN	! if periodic in this dimension
	     IF (x(jdim,i).GT.xmax(jdim)) THEN
!	        print*,' xold,xmax,xnew = ',i,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
	        x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
	     ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!	        print*,' xold,xmin,xnew = ',i,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)	     
                x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
             ENDIF	  
	  ENDIF
       ENDDO
    ENDIF	 

 ENDDO

!
!--set ghost particles if ghost boundaries are used
!	 
 IF (ANY(ibound.GE.2)) CALL set_ghost_particles
!
!--call link list to find neighbours
!
 CALL set_linklist
!
!--calculate density by direct summation
!
 IF (icty.LE.0) CALL iterate_density
!
!--calculate primitive variables from conservative variables
!   
 CALL conservative2primitive
!
!--calculate forces/rates of change using predicted quantities
!	 
 CALL get_rates
!
!--Mid-point Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       pmom(:,i) = pmomin(:,i)
       rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)
       IF (iener.NE.0) en(i) = enin(i)
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i)) 
       alpha(:,i) = alphain(:,i)
       hh(i) = hhin(i)
    ELSE
       x(:,i) = xin(:,i) + dt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       pmom(:,i) = pmomin(:,i) + dt*force(:,i)
       IF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + dt*dhdt(i)
	  IF (hh(i).LE.0.) THEN
	     WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
	     CALL quit
	  ENDIF
       ENDIF
       IF (icty.GE.1) THEN
          rho(i) = rhoin(i) + dt*drhodt(i)
          IF (rho(i).LE.0.) THEN
             WRITE(iprint,*) 'step: rho -ve ',i,rho(i)
    	     CALL quit
          ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + dt*dendt(i)
       IF (iavlim.NE.0) alpha(:,i) = alphain(:,i) + dt*daldt(:,i)	   
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldt(:,i)	  
    ENDIF
 ENDDO	 	
!
!--update density using a full summation every so often
!	 
 IF (MOD(nsteps,ndirect).EQ.0) THEN
    CALL iterate_density
    DO i=1,npart
       rhoin(i) = rho(i)
    ENDDO
 ENDIF
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
!
 IF (ANY(ibound.NE.0)) CALL boundary	! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 dt = min(C_force*dtforce,C_cour*dtcourant)
!
!--need to calculate velocity from particle momentum for next predictor step
!  (no need to call this from output if called here)
! 
 CALL conservative2primitive
  
 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
