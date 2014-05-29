!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a RKF45 integrator
!! At the moment there is no XSPH and no direct summation replacements
!!--------------------------------------------------------------------

!-------------------------------------------------------------------
! Runge-Kutta-Fehlberg Integrator
!-------------------------------------------------------------------
module rkf
  implicit none

!--integrator tolerance
  real, parameter :: rkftol = 1.0E-06

!--5th order RK coefficients
  real, parameter :: b15 = 16./135.
  real, parameter :: b25 = 0
  real, parameter :: b35 = 6656./12825.
  real, parameter :: b45 = 28561./56430.
  real, parameter :: b55 = -9./50.
  real, parameter :: b65 = 2./55.

!--4th order RK coefficients
  real, parameter :: b14 = 25./216.
  real, parameter :: b24 = 0
  real, parameter :: b34 = 1408./2565.
  real, parameter :: b44 = 2197./4104.
  real, parameter :: b54 = -1./5.

  real, parameter :: a21 =  0.25
  real, parameter :: a31 = 3./32.
  real, parameter :: a32 = 9./32.
  real, parameter :: a41 = 1932./2197.
  real, parameter :: a42 = -7200./2197.
  real, parameter :: a43 = 7296./2197.
  real, parameter :: a51 = 439./216.
  real, parameter :: a52 = -8.
  real, parameter :: a53 = 3680./513.
  real, parameter :: a54 = -845./4104.
  real, parameter :: a61 = -8./27.
  real, parameter :: a62 = 2.
  real, parameter :: a63 = -3544./2565.
  real, parameter :: a64 = 1859./4104.
  real, parameter :: a65 = -11./40.

end module rkf


subroutine step ()
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
 integer :: i,ii
 real, dimension(ndimV,npart) :: forcein,dBevoldtin
 real, dimension(npart) :: drhodtin,dhdtin,dendtin,uuin,dpsidtin
 real, dimension(3,npart) :: daldtin
 real :: hdt, dtin
 real, dimension(ndim)  :: xtemp
 real, dimension(ndimV) :: vtemp, Btemp
 real :: dtmin, dtmax, dtnext
 real :: rkferr, rkferrratio, dtfactor, errmhd
 logical :: loop
 real, dimension(ndimV,npart) :: v1,v2,v3,v4,v5,v6
 real, dimension(ndimV,npart) :: f1,f2,f3,f4,f5,f6
 real, dimension(ndimV,npart) :: dB1,dB2,dB3,dB4,dB5,dB6
 real, dimension(npart) :: dh1,dh2,dh3,dh4,dh5,dh6
 real, dimension(npart) :: den1,den2,den3,den4,den5,den6
 real, dimension(npart) :: drho1,drho2,drho3,drho4,drho5,drho6

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
! print *, 'DT = ', dt
 dtin = dt
 dtnext = huge(dtnext)

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

 if (any(iavlim.ne.0)) stop 'rkf not implemented with alpha switches'
 if (idivBzero.ge.2) stop 'rkf not implemented with div B correction'
 if (ihvar.eq.1) stop 'ihvar=1 not implemented with rkf'

!--Loop until RKF timestep gives convergence

 DO WHILE (loop)

!
!RKF step 1
!
 ! xin
 ! velin
 ! forcein
 do i=1,npart
    v1(:,i)  = velin(:,i)
    f1(:,i)  = forcein(:,i)
    dB1(:,i) = dBevoldt(:,i)
    dh1(i)   = dhdtin(i)
    den1(i)  = dendt(i)
    drho1(i) = drhodt(i)
 enddo
!
!--RKF step 2
!
 do i=1,npart
    x(:,i) = xin(:,i) + a21*dt*velin(1:ndim,i)
    vel(:,i) = velin(:,i) + a21*dt*forcein(:,i)

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + a21*dt*dBevoldtin(:,i)
    if (icty.ge.1) rho(i) = rhoin(i) + a21*dt*drhodtin(i)
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + a21*dt*dh1(i)
    endif
    if (iener.ne.0) en(i) = enin(i) + a21*dt*dendtin(i)
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + a21*dt*daldtin(:,i),1.0)
!    if (idivBzero.ge.2) psi(i) = psiin(i) + a21*dt*dpsidtin(i)
 enddo
!
!--calculate all derivatives
!
 call derivs
 do i=1,npart
    v2(:,i)  = vel(:,i)
    f2(:,i)  = force(:,i)
    dB2(:,i) = dBevoldt(:,i)
    dh2(i)   = dhdt(i)
    den2(i)  = dendt(i)
    drho2(i) = drhodt(i)
 enddo
!
!--RK step 3
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a31*v1(1:ndim,i) + a32*v2(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a31*f1(:,i) + a32*f2(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + (a31*dB1(:,i) + a32*dB2(:,i))*dt
    if (icty.ge.1) rho(i) = rhoin(i) + (a31*drho1(i) + a32*drho2(i))*dt
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + (a31*dh1(i) + a32*dh2(i))*dt
    endif
    if (iener.ne.0) en(i) = enin(i) + (a31*den1(i) + a32*den2(i))*dt
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
!    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 do i=1,npart
    v3(:,i)  = vel(:,i)
    f3(:,i)  = force(:,i)
    dB3(:,i) = dBevoldt(:,i)
    dh3(i)   = dhdt(i)
    den3(i)  = dendt(i)
    drho3(i) = drhodt(i)
 enddo
!
!--RK step 4
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a41*v1(1:ndim,i) + a42*v2(1:ndim,i) + a43*v3(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a41*f1(:,i) + a42*f2(:,i) + a43*f3(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + (a41*dB1(:,i) + a42*dB2(:,i) + a43*dB3(:,i))*dt
    if (icty.ge.1) rho(i) = rhoin(i) + (a41*drho1(i) + a42*drho2(i) + a43*drho3(i))*dt
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(a41*dh1(i) + a42*dh2(i) + a43*dh3(i))
    endif
    if (iener.ne.0) en(i) = enin(i) + (a41*den1(i) + a42*den2(i) + a43*den3(i))*dt
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
!    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 do i=1,npart
    v4(:,i)  = vel(:,i)
    f4(:,i)  = force(:,i)
    dB4(:,i) = dBevoldt(:,i)
    dh4(i)   = dhdt(i)
    den4(i)  = dendt(i)
    drho4(i) = drhodt(i)
 enddo
!
!--RK step 5
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a51*v1(1:ndim,i) + a52*v2(1:ndim,i) + a53*v3(1:ndim,i) + a54*v4(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a51*f1(:,i) + a52*f2(:,i) + a53*f3(:,i) + a54*f4(:,i))*dt

    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + (a51*dB1(:,i) + a52*dB2(:,i) + a53*dB3(:,i) + a54*dB4(:,i))*dt
    if (icty.ge.1) rho(i) = rhoin(i) + (a51*drho1(i) + a52*drho2(i) + a53*drho3(i) + a54*drho4(i))*dt
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + (a51*dh1(i) + a52*dh2(i) + a53*dh3(i) + a54*dh4(i))*dt
    endif
    if (iener.ne.0) en(i) = enin(i) + (a51*den1(i) + a52*den2(i) + a53*den3(i) + a54*den4(i))*dt
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
!    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 do i=1,npart
    v5(:,i)  = vel(:,i)
    f5(:,i)  = force(:,i)
    dB5(:,i) = dBevoldt(:,i)
    dh5(i)   = dhdt(i)
    den5(i)  = dendt(i)
    drho5(i) = drhodt(i)
 enddo
!
!--RK step 6
!
 do i=1,npart
    x(:,i) = xin(:,i) + (a61*v1(1:ndim,i) + a62*v2(1:ndim,i) + a63*v3(1:ndim,i) &
         + a64*v4(1:ndim,i) + a65*v5(1:ndim,i))*dt
    vel(:,i) = velin(:,i) + (a61*f1(:,i) + a62*f2(:,i) + a63*f3(:,i) + &
         a64*f4(:,i) + a65*f5(:,i))*dt
    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*(a61*dB1(:,i) + a62*dB2(:,i) + a63*dB3(:,i) + &
         a64*dB4(:,i) + a65*dB5(:,i))*dt
    if (icty.ge.1) rho(i) = rhoin(i) + (a61*drho1(i) + a62*drho2(i) + a63*drho3(i) + a64*drho4(i) + a65*drho5(i))*dt
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + (a61*dh1(i) + a62*dh2(i) + a63*dh3(i) + a64*dh4(i) + a65*dh5(i))*dt
    endif
    if (iener.ne.0) en(i) = enin(i) + (a61*den1(i) + a62*den2(i) + a63*den3(i) + a64*den4(i) + a65*den5(i))*dt
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
!    if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
 enddo
!
!--calculate all derivatives
!
 call derivs
 do i=1,npart
    v6(:,i)  = vel(:,i)
    f6(:,i)  = force(:,i)
    dB6(:,i) = dBevoldt(:,i)
    dh6(i)   = dhdt(i)
    den6(i)  = dendt(i)
    drho6(i) = drhodt(i)
 enddo
!
!--RK final values
!
 rkferr = 0.
 rkferrratio = 1.
 do i=1,npart
!--5th order result
    xtemp(:) = xin(:,i) + dt*(b15*v1(:,i) + b25*v2(:,i) + b35*v3(:,i) + &
         b45*v4(:,i) + b55*v5(:,i) + b65*v6(:,i))
    vtemp(:) = velin(:,i) + dt*(b15*f1(:,i) + b25*f2(:,i) + b35*f3(:,i) + &
         b45*f4(:,i) + b55*f5(:,i) + b65*f6(:,i))
    if (imhd.ne.0 .and. iresist.ne.2) Btemp(:) = Bevolin(:,i) + dt*(b15*dB1(:,i) + b25*dB2(:,i) + b35*dB3(:,i) + &
         b45*dB4(:,i) + b55*dB5(:,i) + b65*dB6(:,i))
!--4th order result
    x(:,i) = xin(:,i) + dt*(b14*v1(:,i) + b24*v2(:,i) + b34*v3(:,i) + &
         b44*v4(:,i) + b54*v5(:,i))
    vel(:,i) = velin(:,i) + dt*(b14*f1(:,i) + b24*f2(:,i) + b34*f3(:,i) + &
         b44*f4(:,i) + b54*f5(:,i))
    if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*(b14*dB1(:,i) + b24*dB2(:,i) + b34*dB3(:,i) + &
         b44*dB4(:,i) + b54*dB5(:,i))
    if (icty.ge.1) rho(i) = rhoin(i) + dt*(b14*drho1(i) + b24*drho2(i) + b34*drho3(i) + &
         b44*drho4(i) + b54*drho5(i))
    if (iener.ne.0) en(i) = enin(i) + dt*(b14*den1(i) + b24*den2(i) + b34*den3(i) + &
         b44*den4(i) + b54*den5(i))
!    print *, 'Vs ', vel(3,i), vtemp(3)
!    rkferr = abs(vtemp(3) - vel(3,i))
!    rkferr = abs(vtemp(3)/vel(3,i) - 1.)
    do ii = 1, ndim
!       rkferr = MAX(rkferr, abs(vtemp(ii) - vel(ii,i)), abs(xtemp(ii) - x(ii,i)))
       rkferr = MAX(rkferr, abs(vtemp(ii)/vel(ii,i) - 1.), abs(xtemp(ii)/x(ii,i)) - 1.)
       if (imhd.ne.0 .and. iresist.ne.2) then
          errmhd = abs(Btemp(ii)/Bevol(ii,i)) - 1.
          !if (errmhd > rkferr) then 
            ! rkferr = errmhd
            ! print*, 'setting rkferr = ',rkferr
          !endif
       endif
    enddo
!    print *, 'rkerr ', rkferr
    rkferrratio = rkferr/rkftol

    if (ihvar.eq.2 .or. ihvar.eq.3) then
       hh(i) = hhin(i) + dt*(b14*dh1(i) + b24*dh2(i) + b34*dh3(i) + b44*dh4(i) + b54*dh5(i))
    endif
!    if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + hdt*(daldt(:,i)+daldtin(:,i)),1.0)
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
     dtnext = min(dt,dtnext)
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
     dt = min(dt*dtfactor,dtnext)
     print *, 'Courant dt ', dt, dtfactor, loop
  endif

!--Loop until RKF timestep gives convergence
  ENDDO

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
