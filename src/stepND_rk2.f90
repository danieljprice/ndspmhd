!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a 2nd order Runge-Kutta algorithm (i.e. midpoint)
!! At the moment there is no XSPH and no direct summation replacements
!!--------------------------------------------------------------------
         
subroutine step
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
 use utils,    only:cross_product3D
!
!--define local variables
!
 implicit none
 integer :: i,j
 real :: hdt,dtstop
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine step'
!
!--set initial quantities
!
 hdt = 0.5*dt
 do i=1,npart
    xin(:,i) = x(:,i)
    velin(:,i) = vel(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
    psiin(i) = psi(i)
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       dustevolin(i) = dustevol(i)
       if (idust.eq.1) deltavin(:,i) = deltav(:,i)
    endif
 enddo
!
!--move everything to the half step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then ! fixed particles
       if (itype(i).eq.itypebnd) then
          j = ireal(i)
          if (j > 0) then
             x(:,i) = xin(:,i) + hdt*velin(1:ndim,j)
          else
             x(:,i) = xin(:,i)
          endif
       endif
       if (imhd.gt.0) Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(i) = dustevolin(i)
          if (idust.eq.1) deltav(:,i) = deltavin(:,i)
       endif
    else
       x(:,i)   = xin(:,i) + hdt*vel(1:ndim,i)
       vel(:,i) = velin(:,i) + hdt*force(:,i)
       if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
       if (icty.ge.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       if (ihvar.eq.1) then
!           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's
       elseif (ihvar.eq.2 .or. ihvar.eq.3) then
          hh(i) = hhin(i) + hdt*dhdt(i)
       endif
       if (iener.ne.0) en(i) = enin(i) + hdt*dendt(i)
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + hdt*daldt(:,i),1.0)
       if (idivBzero.ge.2) psi(i) = psiin(i) + hdt*dpsidt(i) 
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(i) = dustevolin(i) + hdt*ddustevoldt(i)
          if (idust.eq.1) then
             deltav(:,i) = deltavin(:,i) + hdt*ddeltavdt(:,i)
!             if (dustevol(i).gt.0.) then
!                dtstop = Kdrag/(rho(i)*dustevol(i)*(1. - dustevol(i)))
!                deltav(:,i) = deltav(:,i)*exp(-hdt*dtstop)
!             endif
          endif
       endif
    endif
 enddo

!
!--calculate all derivatives at the half step
!
 call derivs
 if (any(ibound.ne.0)) call boundary        ! inflow/outflow/periodic boundary conditions
!
!--Now do the corrector (full) step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then
       if (itype(i).eq.itypebnd) vel(:,i) = velin(:,i)
       if (imhd.gt.0) Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(i) = dustevolin(i)
          if (idust.eq.1) deltav(:,i)  = deltavin(:,i)
       endif
    else
       x(:,i) = xin(:,i) + dt*vel(:,i)
       vel(:,i) = velin(:,i) + dt*force(:,i)
       if (imhd.ne.0) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldt(:,i)
       if (icty.ge.1) rho(i)     = rhoin(i) + dt*drhodt(i)
       if (ihvar.eq.2) then
          hh(i) = hhin(i) + dt*dhdt(i)
          if (hh(i).le.0.) then
             write(iprint,*) 'step: hh -ve ',i,hh(i)
             call quit
          endif
       elseif (ihvar.eq.3) then
          hh(i) = hh(i) + hdt*dhdt(i)
       endif
       if (iener.ne.0)       en(i)      = enin(i) + dt*dendt(i)
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldt(:,i),1.0)
       if (idivbzero.ge.2)   psi(i)     = psiin(i) + dt*dpsidt(i)
       
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(i) = dustevolin(i)  + dt*ddustevoldt(i)
          if (idust.eq.1) then
             deltav(:,i)  = deltavin(:,i)   + dt*ddeltavdt(:,i)
!             if (dustevol(i).gt.0.) then
!                dustfrac(i) = dustevol(i)
!                dtstop = Kdrag/(rho(i)*dustevol(i)*(1. - dustfrac(i)))
!                deltav(:,i) = deltav(:,i)*exp(-dt*dtstop)
!             endif
          endif
       endif
    endif
 enddo
!
!--calculate all derivatives for next step
!
 call derivs

 if (any(ibound.ne.0)) call boundary        ! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 if (.not.dtfixed) dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag,C_force*dtvisc)

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
