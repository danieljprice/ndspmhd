!!--------------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses the time transformed leapfrog algorithm.
!! (actually this is more like the logarithmic Hamiltonian method)
!!
!! NB: It is time-reversible if the density is
!! calculated by direct summation and the equation of state is isothermal or
!! polytropic and obviously also only if there is no dissipation.
!!
!!--------------------------------------------------------------------------------
         
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
 INTEGER :: i
 REAL :: hdt,w1,whalf,dwdt,Omega,epoti,ds0 ! w0 is saved
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
 ds0 = dt0/Omega0     
 hdt = 0.5*ds0/w0 ! dt0 (a constant) and dtscale (1/W) are saved in the module timestep
       
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
       IF (ihvar.EQ.2 .OR. ihvar.EQ.3) hh(i) = hhin(i) + hdt*dhdt(i)
       IF (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)           
    ENDIF
 ENDDO
!
!--calculate all derivatives
!
 call derivs
!
!--get total potential
!
 Omega = 0.
 DO i=1,npart
    call external_potentials(iexternal_force,x(:,i),epoti,ndim)
    Omega = Omega + pmass(i)*(uu(i) + epoti)
 ENDDO
 whalf = 1./Omega
!
!--Corrector step for v (force transformed using omega)
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       vel(:,i) = velin(:,i)
    ELSE
       !--this overwrites the predictor step for the velocity
       vel(:,i) = (velin(:,i) + ds0*force(:,i)/whalf)/(1.+damp) ! stepped through whole timestep
    ENDIF              
 ENDDO 
!
!--get change in time
!
 dwdt = 0.
 DO i=1,npart
    dwdt = dwdt - 0.5*pmass(i)*dot_product(velin(:,i) + vel(:,i),force(:,i))
 ENDDO
 print*,'dwdt = ',dwdt,dwdt*ds0/Omega**2
 dwdt = -dwdt*ds0/Omega**2
!
!--evolve time transformation factor forwards
!
 w1 = w0 + dwdt
!
!--determine new timstep
!
 dtscale = Omega0/w1
 hdt = 0.5*ds0/w1
 
 dt = 0.5*ds0*(1./w0 + 1./w1) ! this is the total time advanced over the whole timestep
 
 print*,'dt 0 ',dt0,'w1 = ',w1,'dtscale = ',dtscale,' dt = ',dt, &
        ' dt/dtc =',dt/(C_cour*dtcourant),' dt(vel) = ',ds0/whalf
!! dthalf = min(C_force*dtforce,C_cour*dtcourant)
 read*
!
!--Corrector step for position
!
 DO i=1,npart
    x(1:ndim,i) = x(1:ndim,i) + hdt*vel(1:ndim,i) ! nb x is stepped from its current value at t^1/2
 ENDDO
!
!--cross boundaries if necessary
!
 IF (ANY(ibound.NE.0)) CALL boundary  ! inflow/outflow/periodic boundary conditions
!
!--calculate drho/dt, den/dt and du/dt using v and x at full timestep
!  and using h, rho, P, B, at half timestep (ie. from previous calculation/call to conservative2primitive)
!
! IF (ANY(ibound.GE.2)) CALL set_ghost_particles
! CALL set_linklist
! CALL ratesvel(x,vel,rho,pr,Bevol)

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
 w0 = w1

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step
