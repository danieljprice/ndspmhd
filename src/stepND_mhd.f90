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
!--Mid-point Predictor step
!      
 hdt = 0.5*dt
! DO i=1,npart
!    xin(:,i) = x(:,i)
!    velin(:,i) = vel(:,i)
!    Bconsin(:,i) = Bcons(:,i)
!    rhoin(i) = rho(i)
!    hhin(i) = hh(i)
!    enin(i) = en(i)
!    alphain(i) = alpha(i)
!    psiin(i) = psi(i)
! ENDDO         
 
 DO i=1,npart
    IF (itype(i).EQ.1) THEN	! fixed particles
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       Bcons(:,i) = Bconsin(:,i)	     
       IF (iener.NE.0) en(i) = enin(i)
       hh(i) = hhin(i)	    
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(i) = alphain(i)
       psi(i) = psiin(i)
    ELSE
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + hdt*dBconsdt(:,i)
       IF (ihvar.EQ.1) THEN
!	   hh(i) = hfact(pmass(i)/rho(i))**hpower	! my version
	  hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower		! Joe's	   
       ELSEIF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (iavlim.NE.0) alpha(i) = alphain(i) + hdt*daldt(i)	 
       psi(i) = psiin(i) + hdt*dpsidt(i)  
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
 CALL link
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
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       Bcons(:,i) = Bconsin(:,i)
       IF (iener.NE.0) en(i) = enin(i)
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i)) 
       alpha(i) = alphain(i)
       hh(i) = hhin(i)
       psi(i) = psiin(i)
    ELSE
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       vel(:,i) = velin(:,i) + hdt*force(:,i)	    
       IF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
	  IF (hh(i).LE.0.) THEN
	     WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
	     CALL quit
	  ENDIF
       ENDIF
       IF (icty.GE.1) THEN
          rho(i) = rhoin(i) + hdt*drhodt(i)
          IF (rho(i).LE.0.) THEN
             WRITE(iprint,*) 'step: rho -ve ',i,rho(i)
    	     CALL quit
          ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (iavlim.NE.0) alpha(i) = alphain(i) + hdt*daldt(i)	   
       IF (imhd.NE.0) Bcons(:,i) = Bconsin(:,i) + hdt*dBconsdt(:,i)
       psi(i) = psiin(i) + hdt*dpsidt(i)	  
    ENDIF 
	      
 ENDDO
	 	
!
!--update density using a full summation every so often
!	 
 IF (MOD(nsteps,ndirect).EQ.0) THEN
    DO i=1,npart
       xin(:,i) = 2.*x(:,i) - xin(:,i)
       x(:,i) = xin(:,i)
    ENDDO
	   
    CALL density
!    ikernavprev = ikernav
!    ikernav = 3
!    CALL iterate_density ! renormalise density *and* smoothing length
!    ikernav = ikernavprev
    IF (ANY(ibound.GT.1)) THEN
       DO i=npart+1,ntotal		! update ghosts
          j = ireal(i)
          rho(i) = rho(j)
          hh(i) = hh(j)
          gradh(i) = gradh(j)       
       ENDDO
    ELSEIF (ANY(ibound.EQ.1)) THEN           ! rewrite over fixed particles
       WHERE (itype(:) .EQ. 1)
          rho(:) = rhoin(:)
	  hh(:) = hhin(:)
       END WHERE	  
    ENDIF
	   
    DO i=1,npart
       rhoin(i) = rho(i)
       velin(:,i) = 2.*vel(:,i)-velin(:,i)
       vel(:,i) = velin(:,i)
       IF (imhd.NE.0) THEN
          Bconsin(:,i) = 2.*Bcons(:,i) - Bconsin(:,i)
          Bcons(:,i) = Bconsin(:,i)
       ENDIF	     
       IF (iener.NE.0) THEN
	  enin(i) = 2.*en(i) - enin(i)
	  en(i) = enin(i)
       ENDIF
       IF (iavlim.NE.0) THEN
	  alphain(i) = 2.*alpha(i) - alphain(i)
          alpha(i) = alphain(i)
       ENDIF
       IF (ihvar.EQ.2) THEN
	  hhin(i) = 2.*hh(i) - hhin(i)
   	  hh(i) = hhin(i)
!         hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower
!	  hhin(i) = hh(i)
       ENDIF	     
       hhin(i) = hh(i)
       psiin(i) = 2.*psi(i) - psiin(i)
       psi(i) = psiin(i)
    ENDDO
	 
 ELSE
	   
    DO i=1,npart
       xin(:,i) = 2.*x(:,i) - xin(:,i)
       x(:,i) = xin(:,i)
       velin(:,i) = 2.*vel(:,i) - velin(:,i)
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) THEN		
          rhoin(i) = 2.*rho(i) - rhoin(i)
	  rho(i) = rhoin(i)
       ELSE
	  rhoin(i) = rho(i) 	 
       ENDIF
       IF (ihvar.EQ.1) THEN		! Joe's version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower
       ELSEIF (ihvar.EQ.2) THEN
          hhin(i) = 2.*hh(i) - hhin(i)
          hh(i) = hhin(i)
          IF (hh(i).LE.0.) THEN
	     WRITE(iprint,*) 'step: corrector: hh -ve ',i,hh(i)
	     CALL quit
          ENDIF       
       ENDIF
       hhin(i) = hh(i)	        
       IF (iavlim.NE.0) THEN
	  alphain(i) = 2.*alpha(i) - alphain(i)
	  alpha(i) = alphain(i)
       ENDIF	     
       IF (imhd.NE.0) THEN
          Bconsin(:,i) = 2.*Bcons(:,i) - Bconsin(:,i)
          Bcons(:,i) = Bconsin(:,i)
       ENDIF	     
       IF (iener.NE.0) THEN
          enin(i) = 2.*en(i) - enin(i)
          en(i) = enin(i)
       ENDIF
       psiin(i) = 2.*psi(i) - psiin(i)
       psi(i) = psiin(i)
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
 
 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
