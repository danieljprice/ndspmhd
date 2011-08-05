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
 use particlesplit, only:particle_splitting
 use geometry, only:coord_transform,vector_transform
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,jdim,ikernavprev,ierr,nerror,nsplit
 REAL, DIMENSION(ndimV,npart) :: forcein,dBevoldtin
 REAL, DIMENSION(npart) :: drhodtin,dhdtin,dendtin,uuin,dpsidtin
 REAL, DIMENSION(3,npart) :: daldtin
 REAL :: hdt
 real, dimension(ndim) :: xcyl,velcyl
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
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    uuin(i) = uu(i)
    alphain(:,i) = alpha(:,i)
    psiin(i) = psi(i)
    
    forcein(:,i) = force(:,i)
    dBevoldtin(:,i) = dBevoldt(:,i)
    drhodtin(i) = drhodt(i)
    dhdtin(i) = dhdt(i)
    dendtin(i) = dendt(i)
    daldtin(:,i) = daldt(:,i)
    dpsidtin(i) = dpsidt(i)
 ENDDO
!
!--if doing divergence correction then do correction to magnetic field
! 
 IF (idivBzero.EQ.10) CALL divBcorrect(npart,ntotal)
 IF (idivBzero.EQ.10) CALL divBcorrect(npart,ntotal)
 IF (idivBzero.EQ.10) CALL divBcorrect(npart,ntotal)

!
!--Leapfrog Predictor step
!      
      
 DO i=1,npart
    IF (itype(i).EQ.1 .or. itype(i).EQ.2) then ! fixed particles
       if (ireal(i).ne.0 .and. itype(i).EQ.1) then
          j = ireal(i)
          x(:,i) = xin(:,i) + dt*velin(1:ndim,j) + 0.5*dt*dt*forcein(1:ndim,j)
       elseif (itype(i).eq.2) then  ! velocities are vr, vphi
          call coord_transform(xin(1:ndim,i),ndim,1,xcyl(:),ndim,2)
          velcyl(1) = 0.
          velcyl(2) = xcyl(1)*omegafixed
          call vector_transform(xcyl(1:ndim),velcyl(:),ndim,2,vel(1:ndim,i),ndim,1)
          x(:,i) = xin(:,i) + dt*vel(:,i)
       else
          write(iprint,*) 'step: error: ireal not set for fixed part ',i,ireal(i)
          stop
       endif
       Bevol(:,i) = Bevolin(:,i)             
       rho(i) = rhoin(i)
       hh(i) = hhin(i)            
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
    ELSE
       x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)           
       vel(:,i) = velin(:,i) + dt*forcein(:,i)
       velin(:,i) = velin(:,i) + 0.5*dt*forcein(:,i)
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i) + dt*drhodtin(i)
       IF (ihvar.EQ.1) THEN
!           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! Joe's           
       ELSEIF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
          hh(i) = hhin(i) + dt*dhdtin(i)
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + dt*dendtin(i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = MIN(alphain(:,i) + dt*daldtin(:,i),1.0)
       IF (idivBzero.GE.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
    ENDIF
 ENDDO
!
!--calculate all derivatives
!
 call derivs
!
!--Leapfrog Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.1 .or. itype(i).EQ.2) THEN
       if (itype(i).EQ.1) vel(:,i) = velin(:,i)
       Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
    ELSE
       vel(:,i) = velin(:,i) + hdt*(force(:,i)) !+forcein(:,i))            
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + hdt*(dBevoldt(:,i)+dBevoldtin(:,i))          
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*(drhodt(i)+drhodtin(i))
       IF (ihvar.EQ.2) THEN
          hh(i) = hhin(i) + hdt*(dhdt(i)+dhdtin(i))
          IF (hh(i).LE.0.) THEN
             WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
             CALL quit
          ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + hdt*(dendt(i)+dendtin(i))
       IF (ANY(iavlim.NE.0)) alpha(:,i) = MIN(alphain(:,i) + hdt*(daldt(:,i)+daldtin(:,i)),1.0)
       IF (idivBzero.GE.2) psi(i) = psiin(i) + hdt*(dpsidt(i)+dpsidtin(i))           
    ENDIF 
              
 ENDDO
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
 IF (ANY(ibound.NE.0)) CALL boundary        ! inflow/outflow/periodic boundary conditions

 IF (isplitpart.GT.0) THEN
    CALL particle_splitting(nsplit)
    IF (nsplit.gt.0) CALL derivs
 ENDIF
!
!--set new timestep from courant/forces condition
!
 dt = min(C_force*dtforce,C_cour*dtcourant)
 !if (C_cour*dtav.lt.dt) then
 !   print*,'WARNING: AV controlling timestep: (old)dt = ',dt,' (new)dt = ',C_cour*dtav
 !   dt = C_cour*dtav
 !endif

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
