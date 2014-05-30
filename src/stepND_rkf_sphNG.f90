!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!
!! This version uses a 2nd order Runge-Kutta-Fehlberg algorithm
!! identical to that used in sphNG
!!
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
 use particlesplit, only:particle_splitting
 use geometry, only:coord_transform,vector_transform
 use utils,    only:cross_product3D
!
!--define local variables
!
 implicit none
 integer :: i,j
 real :: hdt
 real, parameter :: f21 = 1./256., f22 = 255./256.
 real, parameter :: e1 = 1./512.
 real, parameter :: tolpart = 1.e-5, tolhh = 1.e-5
 real :: errx(ndim), errv(ndimV), errh, erru, errB(ndimV), errm, errdivtol, dum2veli(ndimV)
 real :: ratio, dtf21, dtf22, rmod1, rmod, dttol
 real, dimension(ndimV,npart) :: forcein,dBevoldtin,ddeltavdtin
 real, dimension(npart)       :: drhodtin,dhdtin,dendtin,dpsidtin,ddustfracdtin
 real, dimension(3,npart)     :: daldtin
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
       dustfracin(i) = dustfrac(i)
       if (idust.eq.1) deltavin(:,i) = deltav(:,i)
       ddustfracdtin(i) = ddustfracdt(i)
       ddeltavdtin(:,i) = ddeltavdt(:,i)
    endif
    forcein(:,i) = force(:,i)
    drhodtin(i)  = drhodt(i)
    dendtin(i)   = dendt(i)
    dhdtin(i)    = dhdt(i)
    daldtin(:,i) = daldt(:,i)
    dBevoldtin(:,i) = dBevoldt(:,i)
    dpsidtin(i)     = dpsidt(i)
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
          dustfrac(i) = dustfracin(i)
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
          dustfrac(i) = dustfracin(i) + hdt*ddustfracdt(i)
          if (idust.eq.1) then
             deltav(:,i) = deltavin(:,i) + hdt*ddeltavdt(:,i)
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
 dtf21 = dt*f21
 dtf22 = dt*f22
 rmod = huge(rmod)
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
          dustfrac(i) = dustfracin(i)
          if (idust.eq.1) deltav(:,i)  = deltavin(:,i)
       endif
       dum2veli(:) = vel(:,i)
    else
       x(:,i) = xin(:,i) + dtf21*velin(:,i) + dtf22*vel(:,i)
       dum2veli(:) = vel(:,i)
       vel(:,i) = velin(:,i) + dtf21*forcein(:,i) + dtf22*force(:,i)
       if (imhd.ne.0) Bevol(:,i) = Bevolin(:,i) + dtf21*dBevoldtin(:,i) + dtf22*dBevoldt(:,i)
       if (icty.ge.1) rho(i)     = rhoin(i) + dtf21*drhodtin(i) + dtf22*drhodt(i)
       if (ihvar.eq.2) then
          hh(i) = hhin(i) + dtf21*dhdtin(i) + dtf22*dhdt(i)
          if (hh(i).le.0.) then
             write(iprint,*) 'step: hh -ve ',i,hh(i)
             call quit
          endif
       endif
       if (iener.ne.0)       en(i)      = enin(i) + dtf21*dendtin(i) + dtf22*dendt(i)
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dtf21*daldtin(:,i) + dtf22*daldt(:,i),1.0)
       if (idivbzero.ge.2)   psi(i)     = psiin(i) + dtf21*dpsidtin(i) + dtf22*dpsidt(i)
       
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustfrac(i) = dustfracin(i)  + dtf21*ddustfracdtin(i) + dtf22*ddustfracdt(i)
          if (idust.eq.1) then
             deltav(:,i)  = deltavin(:,i)   + dtf21*ddeltavdtin(:,i) + dtf22*ddeltavdt(:,i)
          endif
       endif
    endif
!
!--estimate maximum error
!
    errx(:) = abs(vel(:,i) - dum2veli(:))
    errv(:) = abs(forcein(:,i) - force(:,i))
    erru = abs(dendtin(i) - dendt(i))
    errh = 3.*abs(dhdtin(i) - dhdt(i))
    if (imhd.ne.0) then
       errB(:) = abs(dBevoldtin(:,i) - dBevoldt(:,i))
       errm = max(maxval(errx),maxval(errv),erru,errh,maxval(errB))
    else
       errm = max(maxval(errx),maxval(errv),erru,errh)
    endif
    errdivtol = max(errm/tolpart, errh/tolhh)
!
!--compute ratio of present to next timestep
!
    if (dt > epsilon(dt)) then
       ratio = dt*errdivtol*e1 + epsilon(e1)
       rmod1 = 1./sqrt(ratio)
       !if (rmod1 < rmod) print*,'dt = ',dt,' errm = ',errm, 'errdivtol = ',errdivtol, 'ratio = ',rmod1,rmod
       rmod = min(rmod1,rmod)
    else
       rmod = 1.
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
 if (.not.dtfixed) then
    if (dt > epsilon(dt)) then
       dttol = dt*rmod
    else
       dttol = huge(dttol)
    endif
    print*,' dt cour  = ',C_cour*dtcourant,' dt tol = ',dttol,' old dt = ',dt,' ratio = ',rmod
    dt = min(dttol,C_force*dtforce,C_cour*dtcourant,C_force*dtdrag,C_force*dtvisc)
 endif
 !print*,' setting dt  = ',dt
 !read*

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step