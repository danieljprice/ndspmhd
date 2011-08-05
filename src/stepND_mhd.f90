!!------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!
!! This subroutine is a slightly modified version of the 
!! predictor-corrector scheme. Algorithm goes as follows:
!!
!! Predictor:
!!
!! v^{1/2}   = v^0 + dt/2*f^(-1/2)
!! x^{1/2}   = x^0 + dt/2*v^{1/2}   (note used updated v)
!! rho^{1/2} = rho^0 + dt/2*drhodt^{-1/2}
!!
!!  --> calculate f^{1/2}, drhodt^{1/2} using x^{1/2} and v^{1/2}
!!
!! Corrector:
!!
!! v^*  = v^0 + dt/2*f^{1/2}
!! x^*   = x^0 + dt/2*v^*          (note uses updated v)
!! rho^* = rho^0 + dt/2*drhodt^*
!!
!! v^1   = 2*v^* - v^0     = v^0 + dt*f^{1/2}
!! x^1   = 2*x^* - x^0     = x^0 + dt*v^*
!! rho^1 = 2*rho^* - rho^0 = rho^0 + dt*drhodt^{1/2}
!!
!! Energy, smoothing length, alpha and magnetic field follow density
!!
!! Gives good results (and good stability) on wave and shock-type problems
!! with a very modest courant number (C_cour = 0.8 seems to be enough).
!! On these problems seems to do much better than leapfrog/symplectic,
!! mainly related to the derivatives which depend on velocity (for which the
!! other two methods are not reversible/symplectic).
!!
!! However, results for disks and orbits are slightly better with leapfrog
!!
!!---------------------------------------------------------------------------
         
SUBROUTINE step
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
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
 INTEGER :: i,j,jdim
 REAL :: hdt, maxdivB
 REAL, DIMENSION(ndimV,SIZE(rho)) :: gradpsiprev
 REAL, DIMENSION(SIZE(divB)) :: divBprev
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
!    Bevolin(:,i) = Bevol(:,i)
!    rhoin(i) = rho(i)
!    hhin(i) = hh(i)
!    enin(i) = en(i)
!    alphain(:,i) = alpha(:,i)
!    psiin(i) = psi(i)
! ENDDO         
!
!--if doing divergence correction then do correction to magnetic field
! 
 nsubsteps_divB = 0
 IF (idivBzero.GT.10 .AND. MOD(nsteps,10).EQ.0) THEN
    CALL divBcorrect(npart,ntotal)
    Bevolin(:,1:ntotal) = Bevol(:,1:ntotal)
    IF (any(ibound.ne.0)) WRITE(iprint,*) 'Warning: boundaries not correct'
 ELSEIF (idivBzero.GE.2) THEN
    maxdivB = MAXVAL(ABS(divB(1:npart)))
    nsubsteps_divB = 2
    !IF (maxdivB.gt.0.) nsubsteps_divB = INT(LOG10(maxdivB))
    !print*,'nsubsteps_divB = ',nsubsteps_divB

    !--use grad psi and div B from previous step in substepping
    !  (not ones calculated in this rates call)
    gradpsiprev = gradpsi
    divBprev = divB
 ENDIF
   
 DO i=1,npart
    IF (itype(i).EQ.1) THEN        ! fixed particles
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)             
       IF (iener.NE.0) en(i) = enin(i)
       hh(i) = hhin(i)            
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))

       IF (imhd.NE.0) THEN
          Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
          IF (nsubsteps_divB.EQ.0) THEN
             Bevol(:,i) = Bevol(:,i) + hdt*gradpsi(:,i)
             psi(i) = psiin(i) + hdt*dpsidt(i)
          ENDIF
       ENDIF
       IF (ihvar.EQ.1) THEN
!           hh(i) = hfact(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! Joe's           
       ELSEIF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)
    ENDIF
!
!--for periodic boundaries, allow particles to cross the domain
!  (this is only temporary as it is for the predicted quantity)
!    
    IF (ANY(ibound.EQ.3)) THEN
       DO jdim=1,ndim
          IF (ibound(jdim).EQ.3) THEN        ! if periodic in this dimension
             IF (x(jdim,i).GT.xmax(jdim)) THEN
!                print*,' xold,xmax,xnew = ',i,x(jdim,i),xmax(jdim),xmin(jdim) + x(jdim,i) - xmax(jdim)
                x(jdim,i) = xmin(jdim) + x(jdim,i) - xmax(jdim)
             ELSEIF(x(jdim,i).LT.xmin(jdim)) THEN
!                print*,' xold,xmin,xnew = ',i,x(jdim,i),xmin(jdim),xmax(jdim) + x(jdim,i) - xmin(jdim)             
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
!--do substepping for div B correction
!
 IF (nsubsteps_divB.gt.0) THEN
    print*,' max div B = ',maxval(divB(1:npart)),maxloc(divB(1:npart)),rho(maxloc(divB(1:npart))), &
    ' nsubsteps_divB = ',nsubsteps_divB
    print*,'npart = ',npart, 'ntotal = ',ntotal
    IF (imhd.GE.11) THEN
       CALL substep_divB(dt,nsubsteps_divB,Bevolin(:,1:ntotal), &
            psiin(1:ntotal),divB(1:ntotal),gradpsi(:,1:ntotal), &
            x(:,1:ntotal),hh(1:ntotal),pmass(1:ntotal), &
            itype(1:ntotal),npart,ntotal)
    ELSE
       DO i=1,ntotal
          Bfield(:,i) = Bevolin(:,i)*rho(i)
          gradpsiprev(:,i) = gradpsi(:,i)*rho(i)
       ENDDO
       print*,'max grad psi = ',maxval(gradpsiprev(:,1:ntotal))
       CALL substep_divB(dt,nsubsteps_divB,Bfield(:,1:ntotal), &
            psiin(1:ntotal),divB(1:ntotal),gradpsiprev(:,1:ntotal), &
            x(:,1:ntotal),hh(1:ntotal),pmass(1:ntotal), &
            itype(1:ntotal),npart,ntotal)
       DO i=1,ntotal
          Bevolin(:,i) = Bfield(:,i)/rho(i)
          gradpsi(:,i) = gradpsiprev(:,i)/rho(i)
       ENDDO
    ENDIF
 ENDIF

! print*,'----------------------------'
! do i=17,25
!    print*,i,'force = ',force(:,i)
! enddo
!
!--Mid-point Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)
       IF (iener.NE.0) en(i) = enin(i)
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i)) 
       alpha(:,i) = alphain(:,i)
       hh(i) = hhin(i)
       psi(i) = psiin(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*force(:,i)            
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       IF (ihvar.EQ.2 .OR. (ihvar.EQ.3 .and. itsdensity.eq.0)) THEN
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
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)           
       IF (imhd.NE.0) THEN
          Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
          IF (nsubsteps_divB.eq.0) THEN
             Bevol(:,i) = Bevol(:,i) + hdt*gradpsi(:,i)
             psi(i) = psiin(i) + hdt*dpsidt(i)
          ENDIF          
       ENDIF
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
           
    CALL set_linklist
!    ikernavprev = ikernav
!    ikernav = 3
    CALL iterate_density ! renormalise density *and* smoothing length
!    ikernav = ikernavprev
    IF (ANY(ibound.GT.1)) THEN
       DO i=npart+1,ntotal                ! update ghosts
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
          Bevolin(:,i) = 2.*Bevol(:,i) - Bevolin(:,i)
          Bevol(:,i) = Bevolin(:,i)
       ENDIF             
       IF (iener.NE.0) THEN
          enin(i) = 2.*en(i) - enin(i)
          en(i) = enin(i)
       ENDIF
       IF (ANY(iavlim.NE.0)) THEN
          alphain(:,i) = 2.*alpha(:,i) - alphain(:,i)
          alpha(:,i) = alphain(:,i)
       ENDIF
       IF (ihvar.EQ.2) THEN
          hhin(i) = 2.*hh(i) - hhin(i)
             hh(i) = hhin(i)
!         hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim
!          hhin(i) = hh(i)
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
       IF (ihvar.EQ.1) THEN                ! Joe's version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim
       ELSEIF (ihvar.EQ.2 .or. (ihvar.eq.3 .and. itsdensity.eq.0)) THEN
          hhin(i) = 2.*hh(i) - hhin(i)
          hh(i) = hhin(i)
          IF (hh(i).LE.0.) THEN
             WRITE(iprint,*) 'step: corrector: hh -ve ',i,hh(i)
             CALL quit
          ENDIF       
       ENDIF
       hhin(i) = hh(i)                
       IF (ANY(iavlim.NE.0)) THEN
          alphain(:,i) = 2.*alpha(:,i) - alphain(:,i)
          alpha(:,i) = alphain(:,i)
       ENDIF             
       IF (imhd.NE.0) THEN
          Bevolin(:,i) = 2.*Bevol(:,i) - Bevolin(:,i)
          Bevol(:,i) = Bevolin(:,i)
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
! IF (idivBzero.EQ.1 .and. MOD(nsteps,10.EQ.0) .AND. ALL(ibound.EQ.0)) CALL divBcorrect
!
 IF (ANY(ibound.NE.0)) CALL boundary        ! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 dt = min(C_force*dtforce,C_cour*dtcourant)
 
 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
