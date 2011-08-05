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
 INTEGER :: i,j,jdim,ikernavprev
 REAL :: fhmax, fonh, forcemag,v2i,B2i,frac
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine step'
!
!--Mid-point Predictor step
!      
 B2i = 0.
 v2i = 0.
      
 DO i=1,npart
    IF (itype(i).EQ.1) THEN	! fixed particles
       vel(:,i) = velin(:,i)
       rho(i) = rhoin(i)
       Bfield(:,i) = Bfieldin(:,i)	     
       IF (iener.EQ.3) THEN		! choice of energy equation
          en(i) = enin(i)
          IF (imhd.GE.11) THEN		! B 
              B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
	  ELSEIF (imhd.NE.0) THEN		! B/rho
              B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
	  ENDIF
	  v2i = DOT_PRODUCT(velin(:,i),velin(:,i))		
          uu(i) = en(i) - 0.5*v2i - 0.5*B2i
       ELSE
	  uu(i) = uuin(i)
       ENDIF
       hh(i) = hhin(i)	    
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(i) = alphain(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       IF (imhd.NE.0) Bfield(:,i) = Bfieldin(:,i) + hdt*dBfielddt(:,i)
       IF (ihvar.EQ.1) THEN
!	   hh(i) = hfact(pmass(i)/rho(i))**hpower	! my version
	  hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower		! Joe's	   
       ELSEIF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (iener.EQ.3) THEN
	  en(i) = enin(i) + hdt*dendt(i)
          IF (imhd.GE.11) THEN		! B 
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
	  ELSEIF (imhd.NE.0) THEN		! B/rho
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
	  ENDIF		
	  v2i = DOT_PRODUCT(vel(:,i),vel(:,i))		
	  uu(i) = en(i) - 0.5*v2i - 0.5*B2i
       ELSEIF (iener.GE.1) THEN
  	  uu(i) = uuin(i) + hdt*dudt(i)
       ENDIF
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       IF (iavlim.EQ.1) alpha(i) = alphain(i) + hdt*daldt(i)	   
    ENDIF
!
!--for periodic boundaries, allow particles to cross the domain
!  (this is only temporary as it is for the predicted quantity)
!    
    IF (ibound.EQ.3) THEN
       DO jdim=1,ndim
          IF (x(jdim,i).GT.xmax(jdim)) THEN
!	     print*,' xold,xmax,xnew = ',i,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
	     x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
	  ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!	     print*,' xold,xmin,xnew = ',i,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)	     
             x(jdim,i) = xmax(jdim) - (xmin(jdim) - x(jdim,i))
          ENDIF	  
       ENDDO
    ENDIF	 

 ENDDO

!
!--set ghost particles if ghost boundaries are used
!	 
 IF (ibound.GE.2) CALL set_ghost_particles
!
!--call link list to find neighbours
!
 CALL link
!
!--calculate density by direct summation
!
 IF (icty.LE.0) THEN
    CALL iterate_density
    IF (ibound.GT.1) THEN
       DO i=npart+1,ntotal		! update ghosts
          j = ireal(i)
          rho(i) = rho(j)
          hh(i) = hh(j)
          gradh(i) = gradh(j)       
       ENDDO
    ELSEIF (ibound.EQ.1) THEN           ! rewrite over fixed particles
       WHERE (itype(:) .EQ. 1)
          rho(:) = rhoin(:)
	  hh(:) = hhin(:)
       END WHERE	  
    ENDIF   
 ENDIF   
!
!--set pressure using equation of state
!
 CALL equation_of_state(pr,spsound,uu,rho,gamma,SIZE(rho))
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
       rho(i) = rhoin(i)
       Bfield(:,i) = Bfieldin(:,i)
       IF (iener.EQ.3) THEN		! choice of energy equation
          en(i) = enin(i)
          IF (imhd.GE.11) THEN		! B 
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
          ELSEIF (imhd.NE.0) THEN		! B/rho
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
          ENDIF
          v2i = DOT_PRODUCT(vel(:,i),vel(:,i))	       
	  uu(i) = en(i) - 0.5*v2i - 0.5*B2i
       ELSE
          uu(i) = uuin(i)
       ENDIF
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i)) 
       alpha(i) = alphain(i)
       hh(i) = hhin(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*force(:,i)	    
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
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
       IF (iener.EQ.3) THEN
           en(i) = enin(i) + hdt*dendt(i)
       ELSEIF (iener.GE.1) THEN
           uu(i) = uuin(i) + hdt*dudt(i)
       ENDIF
       IF (iavlim.EQ.1) alpha(i) = alphain(i) + hdt*daldt(i)	   
       IF (imhd.NE.0) Bfield(:,i) = Bfieldin(:,i) + hdt*dBfielddt(:,i)
       
!       IF (i.eq.itemp) THEN
!         PRINT*,'hdt = ',hdt
!	  PRINT*,'f,drhodt,dendt = ',force(:,i),drhodt(i),dendt(i)
!	  PRINT*,'x(',i,') = ',xin(:,i),'+',x(:,i)-xin(:,i)
!	  PRINT*,'vx = ',velin(:,i),'+',vel(:,i)-velin(:,i)
!	  PRINT*,'rho = ',rhoin(i),'+',rho(i)-rhoin(i)
!	  PRINT*,'en = ',enin(i),'+',en(i)-enin(i)
!	  READ*
!       ENDIF	  
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
    IF (ibound.GT.1) THEN
       DO i=npart+1,ntotal		! update ghosts
          j = ireal(i)
          rho(i) = rho(j)
          hh(i) = hh(j)
          gradh(i) = gradh(j)       
       ENDDO
    ELSEIF (ibound.EQ.1) THEN           ! rewrite over fixed particles
       WHERE (itype(:) .EQ. 1)
          rho(:) = rhoin(:)
	  hh(:) = hhin(:)
       END WHERE	  
    ENDIF
	   
    DO i=1,npart
!       hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower
!	hhin(i) = hh(i)
       rhoin(i) = rho(i)
       velin(:,i) = 2.*vel(:,i)-velin(:,i)
       vel(:,i) = velin(:,i)
       IF (imhd.NE.0) THEN
          Bfieldin(:,i) = 2.*Bfield(:,i) - Bfieldin(:,i)
          Bfield(:,i) = Bfieldin(:,i)
       ENDIF	     
       IF (iener.EQ.3) THEN
	  enin(i) = 2.*en(i) - enin(i)
	  en(i) = enin(i)                
          IF (imhd.GE.11) THEN		! B 
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
	  ELSEIF (imhd.NE.0) THEN		! B/rho
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
	  ENDIF
	  v2i = DOT_PRODUCT(vel(:,i),vel(:,i))		
	  uu(i) = en(i) - 0.5*v2i - 0.5*B2i
       ELSEIF (iener.GE.1) THEN
	  uuin(i) = 2.*uu(i) - uuin(i)
	  uu(i) = uuin(i)
       ENDIF
       IF (iavlim.EQ.1) THEN
	  alphain(i) = 2.*alpha(i) - alphain(i)
          alpha(i) = alphain(i)
       ENDIF
       IF (ihvar.EQ.2) THEN
	  hhin(i) = 2.*hh(i) - hhin(i)
   	  hh(i) = hhin(i)
       ENDIF	     
       hhin(i) = hh(i)
    ENDDO
	 
 ELSE
	   
    DO i=1,npart
       xin(:,i) = 2.*x(:,i) - xin(:,i)
       x(:,i) = xin(:,i)
       IF (ihvar.EQ.1) THEN		! Joe's version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**hpower
       ELSEIF (ihvar.EQ.2) THEN
          hhin(i) = 2.*hh(i) - hhin(i)
          hh(i) = hhin(i)
       ENDIF
       hhin(i) = hh(i)	        
       IF (hh(i).LE.0.) THEN
	  WRITE(iprint,*) 'step: corrector: hh -ve ',i,hh(i)
	  CALL quit
       ENDIF


       IF (icty.GE.1) THEN		
           rhoin(i) = 2.*rho(i) - rhoin(i)
	   rho(i) = rhoin(i)
!	   IF (ihvar.NE.0) THEN		! my version
!	      hh(i) = hfact*(pmass(i)/rho(i))**hpower
!	   ENDIF
!	   hhin(i) = hh(i)
       ELSE
	  rhoin(i) = rho(i) 	 
       ENDIF   	    
       velin(:,i) = 2.*vel(:,i) - velin(:,i)
       vel(:,i) = velin(:,i)
       IF (imhd.NE.0) THEN
          Bfieldin(:,i) = 2.*Bfield(:,i) - Bfieldin(:,i)
          Bfield(:,i) = Bfieldin(:,i)
       ENDIF	     
       IF (iener.EQ.3) THEN
          enin(i) = 2.*en(i) - enin(i)
          en(i) = enin(i)
          IF (imhd.GE.11) THEN		! B 
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
	  ELSEIF (imhd.NE.0) THEN		! B/rho
             B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
	  ENDIF
	  v2i = DOT_PRODUCT(vel(:,i),vel(:,i))		
	  uu(i) = en(i) - 0.5*v2i - 0.5*B2i
       ELSE
	  uuin(i) = 2.*uu(i) - uuin(i)
	  uu(i) = uuin(i)
       ENDIF
       IF (iavlim.EQ.1) THEN
	  alphain(i) = 2.*alpha(i) - alphain(i)
	  alpha(i) = alphain(i)
       ENDIF	     
    ENDDO
    
 ENDIF
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
!
!--time step criterion from maximum forces/h
!  dtforce (from step) and dtcourant (from rates) are returned
!  to the main timestepping loop
!
 fhmax = 0.0
 dtforce = 1.e6
 DO i=1,npart
    IF ( ANY(force(:,i).GT.1.e8)) THEN
       WRITE(iprint,*) 'step: force ridiculous ',force(:,i)
       CALL quit
    ENDIF
    forcemag = SQRT(DOT_PRODUCT(force(:,i),force(:,i))) 
    IF (hh(i).LE.0) THEN
       WRITE(iprint,*) 'step: hh(',i,') = ',hh(i)
       CALL quit
    ENDIF
    fonh = forcemag/hh(i)
    IF (fonh.GT.fhmax) fhmax = fonh
 ENDDO

 IF (fhmax.LE.0.) THEN
    WRITE(iprint,*) 'step: fhmax <=0 :',fhmax
    CALL quit
 ENDIF    
 IF (dtforce.GT.0.0) dtforce = SQRT(1./fhmax)
 
 IF (ibound.NE.0) CALL boundary	! inflow/outflow/periodic boundary conditions

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
