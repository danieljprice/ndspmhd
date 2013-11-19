!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a RKF45 integrator
!! At the moment there is no XSPH and no direct summation replacements
!!--------------------------------------------------------------------
         
subroutine step (integratorcheck)
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use eos
 use hterms
 use options
 use part
 use part_in
 use rates
 use timestep
 use setup_params
 use xsph
 use particlesplit, only:particle_splitting
 use geometry, only:coord_transform,vector_transform
 use utils,    only:cross_product3D
 use rkf
!
!--define local variables
!
 implicit none
 integer :: i,j,nsplit, ii
 real, dimension(ndimV,npart) :: forcein,dBevoldtin
 real, dimension(npart) :: drhodtin,dhdtin,dendtin,uuin,dpsidtin
 real, dimension(3,npart) :: daldtin
 real :: hdt, dtin
 real, dimension(ndim)  :: xcyl,velcyl
 real, dimension(ndimV) :: vcrossB
 real, dimension(ndim) :: xtemp
 real, dimension(ndimV) :: vtemp
 real :: dtmin, dtmax
 real :: rkferr, rkferrratio, dtfactor
 logical :: loop
 character (len=*), intent (inout) :: integratorcheck

 if (trim(integratorcheck).eq.'query') then
    integratorcheck = 'rkf'
    return
 endif

 loop = .true.

 dtmin = 1.0E-9
 dtmax = tout
 dtfactor = 1.
 rkferr = 0.
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine step'
!
!--set initial quantities
!
 hdt = 0.5*dt
! print *, 'DT = ', dt
 dtin = dt

 do i=1,npart
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
 enddo
!
!--if doing divergence correction then do correction to magnetic field
! 
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)

!--Loop until RKF timestep gives convergence

 DO WHILE (loop)

!
!RKF step 1
!
 ! xin
 ! velin
 ! forcein
 v1(:,:) = velin(:,:)
 f1(:,:) = forcein(:,:)
 dh1 = dhdtin(:)
!
!--RKF step 2
!
 do i=1,npart
    x(:,i) = xin(:,i) + a21*dt*velin(1:ndim,i)
    vel(:,i) = velin(:,i) + a21*dt*forcein(:,i)

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
    if (ihvar.eq.1) then
       !           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
       hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
    elseif (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*a21*dh1(i)
    endif
    if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 v2(:,:) = vel(:,:)
 f2(:,:) = force(:,:)
 dh2(:) = dhdtin(:)
!
!--RK step 3
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a31*v1(1:ndim,i) + a32*v2(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a31*f1(:,i) + a32*f2(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
    if (ihvar.eq.1) then
       !           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
       hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
    elseif (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(a31*dh1(i) + a32*dh2(i))
    endif
    if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 v3(:,:) = vel(:,:)
 f3(:,:) = force(:,:)
 dh3(:) = dhdtin(:)
!
!--RK step 4
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a41*v1(1:ndim,i) + a42*v2(1:ndim,i) + a43*v3(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a41*f1(:,i) + a42*f2(:,i) + a43*f3(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
    if (ihvar.eq.1) then
       !           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
       hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
    elseif (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(a41*dh1(i) + a42*dh2(i) + a43*dh3(i))
    endif
    if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 v4(:,:) = vel(:,:)
 f4(:,:) = force(:,:)
 dh4(:) = dhdtin(:)
!
!--RK step 5
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a51*v1(1:ndim,i) + a52*v2(1:ndim,i) + a53*v3(1:ndim,i) + a54*v4(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a51*f1(:,i) + a52*f2(:,i) + a53*f3(:,i) + a54*f4(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
    if (ihvar.eq.1) then
       !           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
       hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
    elseif (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(a51*dh1(i) + a52*dh2(i) + a53*dh3(i) + a54*dh4(i))
    endif
    if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 v5(:,:) = vel(:,:)
 f5(:,:) = force(:,:)
 dh5(:) = dhdtin(:)
!
!--RK step 6
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a61*v1(1:ndim,i) + a62*v2(1:ndim,i) + a63*v3(1:ndim,i) &
         + a64*v4(1:ndim,i) + a65*v5(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a61*f1(:,i) + a62*f2(:,i) + a63*f3(:,i) + &
         a64*f4(:,i) + a65*f5(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
    if (ihvar.eq.1) then
       !           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
       hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
    elseif (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(a61*dh1(i) + a62*dh2(i) + a63*dh3(i) + a64*dh4(i) &
            + a65*dh5(i))
    endif
    if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 v6(:,:) = vel(:,:)
 f6(:,:) = force(:,:)
 dh6(:) = dhdtin(:)
!
!--RK final values
!
 rkferr = 0.
 do i=1,npart
!--5th order result
    xtemp(:) = xin(:,i) + dt*(b15*v1(:,i) + b25*v2(:,i) + b35*v3(:,i) + &
         b45*v4(:,i) + b55*v5(:,i) + b65*v6(:,i))
    vtemp(:) = velin(:,i) + dt*(b15*f1(:,i) + b25*f2(:,i) + b35*f3(:,i) + &
         b45*f4(:,i) + b55*f5(:,i) + b65*f6(:,i))
!--4th order result
    x(:,i) = xin(:,i) + dt*(b14*v1(:,i) + b24*v2(:,i) + b34*v3(:,i) + &
         b44*v4(:,i) + b54*v5(:,i))
    vel(:,i) = velin(:,i) + dt*(b14*f1(:,i) + b24*f2(:,i) + b34*f3(:,i) + &
         b44*f4(:,i) + b54*f5(:,i))

!    print *, 'Vs ', vel(3,i), vtemp(3)
!    rkferr = abs(vtemp(3) - vel(3,i))
!    rkferr = abs(vtemp(3)/vel(3,i) - 1.)
    do ii = 1, ndim
!       rkferr = MAX(rkferr, abs(vtemp(ii) - vel(ii,i)), abs(xtemp(ii) - x(ii,i)))
       rkferr = MAX(rkferr, abs(vtemp(ii)/vel(ii,i) - 1.), abs(xtemp(ii)/x(ii,i)) - 1.)
    enddo
!    print *, 'rkerr ', rkferr
    rkferrratio = rkferr/rkftol

    if (icty.ge.1) rho(i) = rhoin(i) + hdt*(drhodt(i)+drhodtin(i))
    if (ihvar.eq.2) then
       hh(i) = hhin(i) + hdt*(dhdt(i)+dhdtin(i))
       if (hh(i).le.0.) then
          write(iprint,*) 'step: hh -ve ',i,hh(i)
          call quit
       endif
    endif
    if (iener.ne.0) en(i) = enin(i) + hdt*(dendt(i)+dendtin(i))
    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + hdt*(daldt(:,i)+daldtin(:,i)),1.0)
  enddo

  if (rkferr.lt.tiny(0.)) then
     dtfactor = 1.0
  else
     dtfactor = (rkftol/(2.*rkferr))**0.25
  endif
!  print *, 'rkferr ', rkferrratio, dtfactor
!  print *, 'nums ', rkferr, rkftol
  if (rkferrratio.gt.1.0) then
     dtfactor = MAX(dtfactor, 0.1)
     print *, 'New dt (down)', dt, dt*dtfactor, dtfactor
     dt = dt*dtfactor
     IF (dt.lt.dtmin) THEN
        print *, 'WARNING: timestep below allowed minimum'
        print *, 'Carrying on anyway with new dt'
     ENDIF
  elseif (dtfactor.gt.1.2) then
     dtfactor = MIN(dtfactor,5.)
     print *, 'New dt ( up )', dt, dt*dtfactor, dtfactor
     dt = dt*dtfactor
     IF (dt.gt.dtmax) THEN
        print *, 'WARNING : dt above allowed maximum. Imposing max dt'
        dt = dtmax
     ENDIF
     loop = .false.
  else
     loop = .false.
  endif

  if (any(ibound.ne.0) .and. .NOT.loop) then
!     print *, 'Setting boundary'
     print *, '------------------------------------------------------------'
     call boundary  !--inflow/outflow/periodic boundary conditions
  endif
!
!--set new timestep from courant/forces condition
!
  if (.NOT.loop.and..NOT.dtfixed) then
     dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag)
     dt = dt*dtfactor
     print *, 'Courant dt ', dt, dtfactor, loop
  endif

!--Loop until RKF timestep gives convergence
  ENDDO

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
