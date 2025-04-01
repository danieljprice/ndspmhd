!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!! This version uses a leapfrog predictor-corrector
!! At the moment there is no XSPH and no direct summation replacements
!! Note that we cannot use leapfrog for the GR code as force.ne.dvel/dt
!!--------------------------------------------------------------------

subroutine step_ind
 use part, only:ibin
 use timestep, only:dt,tout
 implicit none
 integer :: nlevels,level
 real :: dtmax

 nlevels = maxval(ibin)
 level = 1
 dtmax = tout
 !dt = dtmax/2**nlevels
 print*,' calling substep: dt = ',dt,' maxlevel = ',nlevels
 if (abs(dt) <= 0.) then
    call substep_respa(level)
 else
    call recursive_respa(0)
 endif

end subroutine step_ind

recursive subroutine recursive_respa(level)
 use part
 use part_in
 use rates
 use timestep, only:tout,dt
 use bound
 use options
 implicit none
 integer, intent(in) :: level
 integer :: i,j,nlevels
 real :: dti
 !real, allocatable :: dendtin(:)

 nlevels = maxval(ibin)
 dt = 0.0001
 dti = 0.0001/2**level
 print*,'level ',level,' dt = ',dti,' nlevels = ',nlevels

 !dendtin = dendt
 ! push velocity
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2 .or. itype(i).eq.itypebnddust) then
       j = ireal(i)
       vel(:,i) = velin(:,j)
    else
       vel(:,i) = vel(:,i) + 0.5*dti*force_bins(:,i,level+1)
       if (iener.ne.0) en(i) = en(i) + 0.5*dti*dendt(i)
    endif
 enddo

 if (level==nlevels) then
    do i=1,npart
       if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2 .or. itype(i).eq.itypebnddust) then
          j = ireal(i)
          if (j > 0) then
             x(:,i) = x(:,i) + dti*velin(1:ndim,j) !+ 0.5*dt*dt*forcein(1:ndim,j)
          else
             x(:,i) = x(:,i) + dti*velin(:,i)
          endif
       else
          x(:,i) = x(:,i) + dti*vel(:,i)
          hh(i) = hh(i) + dti*dhdt(i)
       endif
    enddo
    call derivs()
    !read*
 else
    call recursive_respa(level+1)
 endif

 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2 .or. itype(i).eq.itypebnddust) then
       j = ireal(i)
       vel(:,i) = velin(:,j)
    else
       vel(:,i) = vel(:,i) + 0.5*dti*force_bins(:,i,level+1)
       if (iener.ne.0) en(i) = en(i) + 0.5*dti*dendt(i)
    endif
 enddo
 

end subroutine recursive_respa

subroutine substep_respa(level)
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
!
!--define local variables
!
 implicit none
 integer, intent(in) :: level
 integer :: i,j
 real, dimension(ndimV,npart) :: forcein,dBevoldtin,ddeltavdtin
 real, dimension(npart)       :: drhodtin,dhdtin,dendtin,dpsidtin
 real, dimension(ndust,npart) :: ddustevoldtin
 real, dimension(3,npart,maxlevels) :: force_bins_in
 real, dimension(3,npart) :: daldtin
 real :: hdt,tol,errv,errB,errmax,errvmax,errBmax,dttol
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine step'
 if (imhd.lt.0) stop 'step routine not compatible with vector potential'
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

    forcein(:,i) = force(:,i)
    force_bins_in(:,i,:) = force_bins(:,i,:)
    dBevoldtin(:,i) = dBevoldt(:,i)
    drhodtin(i) = drhodt(i)
    dhdtin(i) = dhdt(i)
    dendtin(i) = dendt(i)
    daldtin(:,i) = daldt(:,i)
    dpsidtin(i) = dpsidt(i)

    if (onef_dust) then
       select case(idustevol)
       case(1)
          dustevolin(:,i) = sqrt(dustevol(:,i)**2)
       case default
          dustevolin(:,i) = dustevol(:,i)
       end select
       ddustevoldtin(:,i) = ddustevoldt(:,i)
       if (idust.eq.1) then
          deltavin(:,i)     = deltav(:,i)
          ddeltavdtin(:,i)  = ddeltavdt(:,i)
       endif
    endif
 enddo
!
!--Leapfrog Predictor step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2 .or. itype(i).eq.itypebnddust) then ! fixed particles
       if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnddust) then
          j = ireal(i)
          if (j > 0) then
             x(:,i) = xin(:,i) + dt*velin(1:ndim,j) + 0.5*dt*dt*forcein(1:ndim,j)
          else
             x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)
          endif
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
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(:,i) = dustevolin(:,i)
          if (idust.eq.1) deltav(:,i) = deltavin(:,i)
       endif
    else
       vel(:,i) = velin(:,i) + 0.5*dt*force_bins(:,i,level)
       x(:,i) = xin(:,i) + dt*vel(1:ndim,i)
       vel(:,i) = vel(:,i) + 0.5*dt*force_bins(:,i,level)

       if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
       if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
       if (ihvar.eq.1) then
!           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's
       elseif (ihvar.eq.2 .or. ihvar.eq.3) then
          hh(i) = hhin(i) + dt*dhdtin(i)
       endif
       if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
       where(iavlim.ne.0) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
       if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i)
       if (onef_dust) then
          dustevol(:,i) = dustevolin(:,i) + dt*ddustevoldtin(:,i)
          if (idust.eq.1) deltav(:,i) = deltavin(:,i) + dt*ddeltavdtin(:,i)
       endif
    endif
 enddo
!
!--calculate all derivatives
!
 force_bins(:,:,level) = 0.
 call derivs
!
!--Leapfrog Corrector step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2 .or. itype(i).eq.itypebnddust) then
       if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnddust) vel(:,i) = velin(:,i)
       Bevol(:,i) = Bevolin(:,i)
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustevol(:,i) = dustevolin(:,i)
          if (idust.eq.1) deltav(:,i)  = deltavin(:,i)
       endif
    else
       vel(:,i) = vel(:,i) + hdt*(force_bins(:,i,level) - force_bins_in(:,i,level))
       if (imhd.ne.0) then
          if (iresist.eq.2) then
             Bevol(:,i) = Bevolin(:,i) + dt*dBevoldt(:,i)
          else
             Bevol(:,i) = Bevolin(:,i) + hdt*(dBevoldt(:,i)+dBevoldtin(:,i))
          endif
       endif
       if (icty.ge.1) rho(i) = rhoin(i) + hdt*(drhodt(i)+drhodtin(i))
       if (ihvar.eq.2) then
          hh(i) = hhin(i) + hdt*(dhdt(i)+dhdtin(i))
          if (hh(i).le.0.) then
             write(iprint,*) 'step: hh -ve ',i,hh(i)
             call quit
          endif
       endif
       if (iener.ne.0) en(i) = enin(i) + hdt*(dendt(i)+dendtin(i))
       where(iavlim.ne.0) alpha(:,i) = min(alphain(:,i) + hdt*(daldt(:,i)+daldtin(:,i)),1.0)
       if (idivbzero.ge.2) psi(i) = psiin(i) + hdt*(dpsidt(i)+dpsidtin(i))

       if (onef_dust) then
          dustevol(:,i) = dustevolin(:,i) + hdt*(ddustevoldt(:,i) + ddustevoldtin(:,i))
          if (idust.eq.1) deltav(:,i) = deltavin(:,i) + hdt*(ddeltavdt(:,i) + ddeltavdtin(:,i))
       endif
    endif

 enddo

 if (any(ibound.ne.0)) call boundary        ! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 if (.not.dtfixed) then
    !
    !--set a tolerance-based timestep based on the difference between the
    !  second order corrector step and the first order prediction
    !
    tol = 1.e-3
    errmax = 0.
    errvmax = 0.
    errBmax = 0.
    do i=1,npart
       errv = 0.5*dt*maxval(abs(force(:,i) - forcein(:,i)))
       errB = 0.5*dt*maxval(abs(dBevoldt(:,i) - dBevoldtin(:,i)))
       errmax = max(errv,errB,errmax)
       errvmax = max(errv,errvmax)
       errBmax = max(errB,errBmax)
    enddo
    if (errmax > epsilon(errmax) .and. dt > epsilon(dt)) then
       dttol = dt*sqrt(tol/errmax)
    else
       dttol = huge(dttol)
    endif

    dt = min(C_force*dtforce,C_cour*dtcourant,0.9*dtdrag,C_force*dtvisc)
 endif

 if (trace) write (iprint,*) ' Exiting subroutine step'

end subroutine substep_respa
