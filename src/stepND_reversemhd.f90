!!--------------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a reversible (symplectic) verlet predictor-corrector method
!!
!! NB: It is time-reversible if the density is
!! calculated by direct summation and the equation of state is isothermal or
!! polytropic and obviously also only if there is no dissipation.
!!
!! At the moment there is no XSPH and no direct summation replacement of density
!!
!!--------------------------------------------------------------------------------
         
SUBROUTINE step (integratorcheck)
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
 INTEGER :: i
 REAL :: hdt, dt1, dthalf
 character (len=*), intent (inout) :: integratorcheck

 if (trim(integratorcheck).eq.'query') then
    integratorcheck = 'reversemhd'
    return
 endif

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
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
 ENDDO
!
!--Predict to the half time step
!
!  in the following we reserve dt to denote the total time advanced through the timestep
!  -> this is 0.5*(dt_0 + dt_1) and is the amount of time incremented in evolve
!
      
 hdt = 0.5*dt0                ! dt0 is saved from the last timestep in the module timestep
       !! print*,'start: dBevoldt = ',dBevoldt(:,1:10)
 dt = 0.0001
 if (nsteps.gt.1000) dt = -dt
 hdt = 0.5*dt

 DO i=1,npart
    IF (itype(i).EQ.1) THEN        ! fixed particles
       x(1:ndim,i) = xin(1:ndim,i) + hdt*velin(1:ndim,i)       
       vel(:,i) = velin(:,i)
       Bevol(:,i) = Bevolin(:,i)             
       rho(i) = rhoin(i)
       hh(i) = hhin(i)            
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
    ELSE
!
!--x, rho and en are predicted to the half timestep
!
       x(1:ndim,i) = xin(1:ndim,i) + hdt*velin(1:ndim,i)
       !!--vel is only needed at this stage for the viscosity (do not damp)
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (ihvar.EQ.2) hh(i) = hhin(i) + hdt*dhdt(i)
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)           
    ENDIF
 ENDDO
!
!--calculate all derivatives
!
 call derivs
 
 !!CALL output(time,nsteps)
!
!--rates also returns the timestep condition at the half time step
!  we work out the new timestep by setting this timestep condition as an
!  average of the old and new timestep: ie.
!  1/dt_1/2 = 1/2 ( 1/dt_0 + 1/dt_1 ) -> use this to get dt_1 given dt_1/2 and dt_0
!
 dthalf = min(C_force*dtforce,C_cour*dtcourant)
 dt1 = dthalf**2/dt0

 !!dt1 = 1./(2./dthalf - 1./dt0)!!! 2.*dthalf - dt0   !!
 IF (dt1.GT.dtcourant) THEN
    WRITE(iprint,*) 'WARNING: dt1>dtcourant, re-syncing steps ',dt1,dtcourant
    dt1 = dthalf
 ENDIF
 
!! dt = 0.5*(dt0 + dt1) ! this is the total time advanced over the whole timestep
 hdt = 0.5*dt  !!0.5*dt1

! print*,'dt 0 ',dt0,' dt1 = ',dt1, 'dthalf = ',dthalf,' 0.5*(dt0 + dt1) = ',dt

!
!--Corrector step (step from the half timestep to the full timestep using the new dt)
!  so total time advanced is 0.5*dt_0 + 0.5*dt_1
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       vel(:,i) = velin(:,i)
       x(1:ndim,i) = x(1:ndim,i) + hdt*vel(1:ndim,i)
    ELSE
       !--this overwrites the predictor step for the velocity
       vel(:,i) = (velin(:,i) + dt*force(:,i))/(1.+damp) ! stepped through whole timestep
       x(1:ndim,i) = x(1:ndim,i) + hdt*vel(1:ndim,i) ! nb x is stepped from its current value at t^1/2
    ENDIF              
 ENDDO
!
!--cross boundaries if necessary
!
 IF (ANY(ibound.NE.0)) CALL boundary  ! inflow/outflow/periodic boundary conditions
!
!--calculate drho/dt, den/dt and du/dt using v and x at full timestep
!  and using h, rho, P, B, at half timestep (ie. from previous calculation/call to conservative2primitive)
!
 if (imhd.eq.5 .or. iener.eq.2) then
    IF (ANY(ibound.GE.2)) CALL set_ghost_particles
    CALL set_linklist
    DO i=1,npart
       Bfield(:,i) = Bevol(:,i)*rho(i)
    ENDDO
    DO i=npart+1,ntotal
       Bfield(:,i) = Bfield(:,ireal(i))
    ENDDO
    call iterate_density
 endif

 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
    ELSE
       IF (icty.GE.1) rho(i) = rho(i) + hdt*drhodt(i)
       IF (iener.NE.0) en(i) = en(i) + hdt*dendt(i)
       IF (ihvar.EQ.2) hh(i) = hh(i) + hdt*dhdt(i)
       IF (imhd.NE.0) Bevol(:,i) = Bevol(:,i) + hdt*dBevoldt(:,i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)
    ENDIF
 ENDDO
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
