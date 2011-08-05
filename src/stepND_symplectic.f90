!!-------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a reversible (symplectic) verlet predictor-corrector method
!!
!! NB: It is symplectic and manifestly time-reversible if the density is
!! calculated by direct summation and the equation of state is isothermal or
!! polytropic and obviously also only if there is no dissipation.
!!
!! At the moment there is no XSPH and no direct summation replacement of density
!!
!! current version not as good as the leapfrog integrator for anything 
!! involving viscosity
!!-------------------------------------------------------------------------
	 
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
 REAL :: hdt, dt1, dthalf
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine step'
!
!--set initial quantities
!
 DO i=1,npart
    xin(:,i) = x(:,i)
    velin(:,i) = vel(:,i)
    Bconsin(:,i) = Bcons(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(i) = alpha(i)
 ENDDO
!
!--Predict to the half time step
!
!  in the following we reserve dt to denote the total time advanced through the timestep
!  -> this is 0.5*(dt_0 + dt_1) and is the amount of time incremented in evolve
!
      
 hdt = 0.5*dt0		! dt0 is saved from the last timestep in the module timestep
       
 DO i=1,npart
    IF (itype(i).EQ.1) THEN	! fixed particles
       x(:,i) = xin(:,i) + hdt*velin(1:ndim,i)       
       vel(:,i) = velin(:,i)
       Bcons(:,i) = Bconsin(:,i)	     
       rho(i) = rhoin(i)
       hh(i) = hhin(i)	    
       en(i) = enin(i)
       alpha(i) = alphain(i)
    ELSE
!
!--x is predicted to the half timestep
!
       x(:,i) = xin(:,i) + hdt*velin(1:ndim,i)          
!
!--also update the other quantities to the half timestep - this is not symplectic,
!  however only x should be necessary to calculate the forces if the density is
!  calculated by summation, there is no dissipation (and hence no energy equation)
!  and also no mag field. Could perhaps work out a better way to do this.
!
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + hdt*dBconsdt(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (ihvar.EQ.1) THEN
!	   hh(i) = hfact*(pmass(i)/rho(i))**hpower	! my version
	  hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower		! Joe's	   
       ELSEIF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (iavlim.NE.0) alpha(i) = alphain(i) + hdt*daldt(i)	   
    ENDIF
 ENDDO
!
!--allow particles to cross boundaries (this is final as we update x again from the predicted x)
!
 IF (ANY(ibound.NE.0)) CALL boundary	! inflow/outflow/periodic boundary conditions!
!
!--set ghost particles if ghost boundaries are used
!	 
 IF (ANY(ibound.GE.2)) CALL set_ghost_particles
!
!--call link list to find neighbours
!
 CALL link
!
!--calculate density by direct summation
!
 IF (icty.LE.0) CALL iterate_density
!
!--calculate primitive variables (u,B) from conservative variables (en,B/rho)
!   
 CALL conservative2primitive
!
!--calculate forces/rates of change using predicted quantities
!	 
 CALL get_rates
!
!--rates also returns the timestep condition at the half time step
!  we work out the new timestep by setting this timestep condition as an
!  average of the old and new timestep: ie.
!  1/dt_1/2 = 1/2 ( 1/dt_0 + 1/dt_1 ) -> use this to get dt_1 given dt_1/2 and dt_0
!
 dthalf = min(dtforce,dtcourant)
 dt1 = 2.*dthalf - dt0	! this averaging so dt0 can be = 0 to begin with  !1./(2./dthalf - 1./dt0)
 dt = 0.5*(dt1 + dt0)	! this is the total time advanced over the whole timestep
 hdt = 0.5*dt1
!
!--Corrector step (step from the half timestep to the full timestep using the new dt)
!  so total time advanced is 0.5*dt_0 + 0.5*dt_1
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       vel(:,i) = velin(:,i)
       x(:,i) = x(:,i) + hdt*vel(:,i)
       Bcons(:,i) = Bconsin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(i) = alphain(i)
    ELSE
!--this should be the only update for the velocity (ie. no predictor)
       vel(:,i) = velin(:,i) + dt*force(:,i)	! stepped through whole timestep
       x(:,i) = x(:,i) + hdt*vel(:,i)		! nb x is stepped from its current value at t^1/2
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + dt*dBconsdt(:,i)   ! ** CHECK THE REST **
       IF (icty.GE.1) rho(i) = rhoin(i) + dt*drhodt(i)
       IF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + dt*dhdt(i)
	  IF (hh(i).LE.0.) THEN
	     WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
	     CALL quit
	  ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + dt*dendt(i)
       IF (iavlim.NE.0) alpha(i) = alphain(i) + dt*daldt(i)
    ENDIF 
	      
 ENDDO
!
!--cross boundaries if necessary
!
 IF (ANY(ibound.NE.0)) CALL boundary	! inflow/outflow/periodic boundary conditions
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
!
!
!--reset dt0
!
 dt0 = dt1

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
