!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a leapfrog predictor-corrector
!! At the moment there is no XSPH and no direct summation replacements
!! Note that we cannot use leapfrog for the GR code as force.ne.dvel/dt
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
 INTEGER :: i,j,jdim,ikernavprev,ierr,nerror
 REAL, DIMENSION(ndimV,npart) :: forcein,dBconsdtin
 REAL, DIMENSION(npart) :: drhodtin,dhdtin,dendtin,daldtin,uuin
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
    velin(:,i) = vel(:,i)
    Bconsin(:,i) = Bcons(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    uuin(i) = uu(i)
    alphain(i) = alpha(i)
    
    forcein(:,i) = force(:,i)
    dBconsdtin(:,i) = dBconsdt(:,i)
    drhodtin(i) = drhodt(i)
    dhdtin(i) = dhdt(i)
    dendtin(i) = dendt(i)
    daldtin(i) = daldt(i)
 ENDDO
!
!--Leapfrog Predictor step
!      
      
 DO i=1,npart
    IF (itype(i).EQ.1) THEN	! fixed particles
       x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)       
       vel(:,i) = velin(:,i)
       Bcons(:,i) = Bconsin(:,i)	     
       rho(i) = rhoin(i)
       hh(i) = hhin(i)	    
       en(i) = enin(i)
       alpha(i) = alphain(i)
    ELSE
       x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)           
       vel(:,i) = velin(:,i) + dt*forcein(:,i)
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + dt*dBconsdtin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i) + dt*drhodtin(i)
       IF (ihvar.EQ.1) THEN
!	   hh(i) = hfact*(pmass(i)/rho(i))**hpower	! my version
	  hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower		! Joe's	   
       ELSEIF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
          hh(i) = hhin(i) + dt*dhdtin(i)
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + dt*dendtin(i)
       IF (iavlim.NE.0) alpha(i) = alphain(i) + dt*daldtin(i)	   
    ENDIF
 ENDDO
!
!--allow particles to cross boundary
!
 IF (ANY(ibound.NE.0)) CALL boundary	! inflow/outflow/periodic boundary conditions
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
!--Leapfrog Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       vel(:,i) = velin(:,i)
       Bcons(:,i) = Bconsin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(i) = alphain(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*(force(:,i)+forcein(:,i))	    
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + hdt*(dBconsdt(:,i)+dBconsdtin(:,i))	  
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*(drhodt(i)+drhodtin(i))
       IF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*(dhdt(i)+dhdtin(i))
	  IF (hh(i).LE.0.) THEN
	     WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
	     CALL quit
	  ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + hdt*(dendt(i)+dendtin(i))
       IF (iavlim.NE.0) alpha(i) = alphain(i) + hdt*(daldt(i)+daldtin(i))	   
    ENDIF 
	      
 ENDDO
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
!
!
!--set new timestep from courant/forces condition
!
 dt = min(dtforce,dtcourant)

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
