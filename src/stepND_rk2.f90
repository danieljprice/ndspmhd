!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a leapfrog predictor-corrector
!! At the moment there is no XSPH and no direct summation replacements
!! Note that we cannot use leapfrog for the GR code as force.ne.dvel/dt
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
 use particlesplit, only:particle_splitting
 use geometry, only:coord_transform,vector_transform
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
    if (idust.eq.1) then
       dusttogasin(i)    = dusttogas(i)
       deltavin(:,i)     = deltav(:,i)
    endif
 enddo
!
!--move everything to the half step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then ! fixed particles
       if (ireal(i).ne.0 .and. itype(i).eq.itypebnd) then
          j = ireal(i)
          x(:,i) = xin(:,i) + hdt*velin(1:ndim,j)
       else
          write(iprint,*) 'step: error: ireal not set for fixed part ',i,ireal(i)
          stop
       endif
       if (imhd.gt.0) Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)            
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i)
          deltav(:,i) = deltavin(:,i)
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
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i) + hdt*ddusttogasdt(i)
          deltav(:,i) = deltavin(:,i) + hdt*ddeltavdt(:,i)
          if (dusttogas(i).gt.0.) then
             dtstop = Kdrag*(1. + dusttogas(i))**2/(rho(i)*dusttogas(i))
             deltav(:,i) = deltav(:,i)*exp(-hdt*dtstop)
          endif
       endif
    endif
 enddo

 !if (any(ibound.ne.0)) call boundary        ! inflow/outflow/periodic boundary conditions
!
!--calculate all derivatives at the half step
!
 call derivs
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
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i)
          deltav(:,i)  = deltavin(:,i)
       endif
    else
       vel(:,i) = velin(:,i) + dt*force(:,i)
       if (imhd.ne.0) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldt(:,i)
       if (icty.ge.1) rho(i)     = rhoin(i) + dt*drhodt(i)
       if (ihvar.eq.2) then
          hh(i) = hhin(i) + dt*dhdt(i)
          if (hh(i).le.0.) then
             write(iprint,*) 'step: hh -ve ',i,hh(i)
             call quit
          endif
       endif
       if (iener.ne.0)       en(i)      = enin(i) + dt*dendt(i)
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldt(:,i),1.0)
       if (idivbzero.ge.2)   psi(i)     = psiin(i) + dt*dpsidt(i)
       
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i)  + dt*ddusttogasdt(i)
          deltav(:,i)  = deltavin(:,i)   + dt*ddeltavdt(:,i)
          if (dusttogas(i).gt.0.) then
             dtstop = Kdrag*(1. + dusttogas(i))**2/(rho(i)*dusttogas(i))
             deltav(:,i) = deltav(:,i)*exp(-dt*dtstop)
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
 if (.not.dtfixed) dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag)

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
