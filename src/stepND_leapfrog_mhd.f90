!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
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
 integer :: i,j,nsplit
 real, dimension(ndimV,npart) :: forcein,dBevoldtin,ddeltavdtin
 real, dimension(npart) :: drhodtin,dhdtin,dendtin,dpsidtin,ddustfracdtin
 real, dimension(3,npart) :: daldtin
 real :: hdt,tol,errv,errB,errmax,errvmax,errBmax,dttol
 real, dimension(ndim)  :: xcyl,velcyl
 real, dimension(ndimV) :: vcrossB
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
    
    forcein(:,i) = force(:,i)
    dBevoldtin(:,i) = dBevoldt(:,i)
    drhodtin(i) = drhodt(i)
    dhdtin(i) = dhdt(i)
    dendtin(i) = dendt(i)
    daldtin(:,i) = daldt(:,i)
    dpsidtin(i) = dpsidt(i)
    
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       dustfracin(i)    = dustfrac(i)
       ddustfracdtin(i) = ddustfracdt(i)
       if (idust.eq.1) then
          deltavin(:,i)     = deltav(:,i)
          ddeltavdtin(:,i)  = ddeltavdt(:,i)
       endif
    endif
 enddo
!
!--if doing divergence correction then do correction to magnetic field
! 
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)
!
!--Leapfrog Predictor step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then ! fixed particles
       if (itype(i).eq.itypebnd) then
          j = ireal(i)
          if (j > 0) then
             x(:,i) = xin(:,i) + dt*velin(1:ndim,j) + 0.5*dt*dt*forcein(1:ndim,j)
          else
             x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)
          endif
       elseif (itype(i).eq.itypebnd2) then  ! velocities are vr, vphi
          call coord_transform(xin(1:ndim,i),ndim,1,xcyl(:),ndim,2)
          velcyl(1) = 0.
          velcyl(2) = xcyl(1)*omegafixed
          call vector_transform(xcyl(1:ndim),velcyl(:),ndim,2,vel(1:ndim,i),ndim,1)
          x(:,i) = xin(:,i) + dt*vel(:,i)
       else
          write(iprint,*) 'step: error: ireal not set for fixed part ',i,ireal(i)
          stop
       endif
       if (imhd.lt.0) then
          call cross_product3D(velin(:,i),Bconst(:),vcrossB)
          Bevol(:,i) = Bevolin(:,i) + dt*vcrossB(:)
          dBevoldtin(:,i) = vcrossB(:)
       else
          Bevol(:,i) = Bevolin(:,i)
       endif
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
       x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)           
       vel(:,i) = velin(:,i) + dt*forcein(:,i)
       if (imhd.ne.0 .and. iresist.ne.2) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldtin(:,i)
       if (icty.ge.1) rho(i) = rhoin(i) + dt*drhodtin(i)
       if (ihvar.eq.1) then
!           hh(i) = hfact*(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! joe's           
       elseif (ihvar.eq.2 .or. ihvar.eq.3) then
          hh(i) = hhin(i) + dt*dhdtin(i)
       endif
       if (iener.ne.0) en(i) = enin(i) + dt*dendtin(i)
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + dt*daldtin(:,i),1.0)
       if (idivBzero.ge.2) psi(i) = psiin(i) + dt*dpsidtin(i) 
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustfrac(i) = dustfracin(i) + dt*ddustfracdtin(i)
          if (idust.eq.1) deltav(:,i) = deltavin(:,i) + dt*ddeltavdtin(:,i)
       endif
    endif
 enddo
!
!--calculate all derivatives
!
 call derivs
!
!--Leapfrog Corrector step
!
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then
       if (itype(i).eq.itypebnd) vel(:,i) = velin(:,i)
       if (imhd.lt.0) then
          call cross_product3D(vel(:,i),Bconst(:),vcrossB)
          Bevol(:,i) = Bevolin(:,i) + hdt*(vcrossB(:) + dBevoldtin(:,i))
       else
          Bevol(:,i) = Bevolin(:,i)
       endif
       rho(i) = rhoin(i)
       hh(i) = hhin(i)
       en(i) = enin(i)
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustfrac(i) = dustfracin(i)
          if (idust.eq.1) deltav(:,i)  = deltavin(:,i)
       endif
    else
       vel(:,i) = velin(:,i) + hdt*(force(:,i) + forcein(:,i))            
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
       if (any(iavlim.ne.0)) alpha(:,i) = min(alphain(:,i) + hdt*(daldt(:,i)+daldtin(:,i)),1.0)
       if (idivbzero.ge.2) psi(i) = psiin(i) + hdt*(dpsidt(i)+dpsidtin(i))           
       
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustfrac(i) = dustfracin(i) + hdt*(ddustfracdt(i) + ddustfracdtin(i))
          if (idust.eq.1) deltav(:,i) = deltavin(:,i) + hdt*(ddeltavdt(:,i) + ddeltavdtin(:,i))
       endif
    endif

 enddo
!
!--if doing divergence correction then do correction to magnetic field
! 
! IF (idivBzero.NE.0) CALL divBcorrect
 if (any(ibound.ne.0)) call boundary        ! inflow/outflow/periodic boundary conditions

 if (isplitpart.gt.0) then
    call particle_splitting(nsplit)
    if (nsplit.gt.0) call derivs
 endif
!
!--set new timestep from courant/forces condition
!
 if (.not.dtfixed) then
    !
    !--set a tolerance-based timestep based on the difference between the
    !  second order corrector step and the first order prediction
    !
    tol = 1.e-4
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
    
    dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag,C_force*dtvisc)
    if (dttol < dt) then
       dt = dttol
       print "(5(a,es10.3))",'dt (tol) = ',dt,' fac=',sqrt(tol/errmax),' Errmax = ',errmax,' Err v:',errvmax,' Err B:',errBmax
    endif
    !if (C_cour*dtav.lt.dt) then
    !   print*,'WARNING: AV controlling timestep: (old)dt = ',dt,' (new)dt = ',C_cour*dtav
    !   dt = C_cour*dtav
    !endif
    !dtrho = huge(dtrho)
    !do i=1,npart
    !   dtrhoi = abs(rho(i)/(drhodt(i) + epsilon(drhodt)))
    !   dtrho = min(dtrho,0.01*dtrhoi)
    !enddo
    !if (C_rho*dtrho/dtcourant .lt. C_cour) then
    !   dt = dtrho
    !   write(iprint,*) 'dtrho equiv courant number = ',C_rho*dtrho/dtcourant
    !endif
 endif

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
