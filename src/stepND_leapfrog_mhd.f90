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
 real, dimension(npart) :: drhodtin,dhdtin,dendtin,uuin,dpsidtin,ddusttogasdtin
 real, dimension(3,npart) :: daldtin
 real :: hdt
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
    
    if (idust.eq.1) then
       ddusttogasdtin(i) = ddusttogasdt(i)
       ddeltavdtin(:,i)  = ddeltavdt(:,i)
    endif
 enddo
!
!--if doing divergence correction then do correction to magnetic field
! 
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)
 if (idivbzero.eq.10) call divBcorrect(npart,ntotal)

!
!--Leapfrog Predictor step
!      
      
 do i=1,npart
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnd2) then ! fixed particles
       if (ireal(i).ne.0 .and. itype(i).eq.itypebnd) then
          j = ireal(i)
          x(:,i) = xin(:,i) + dt*velin(1:ndim,j) + 0.5*dt*dt*forcein(1:ndim,j)
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
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i)
          deltav(:,i) = deltavin(:,i)
       endif
    else
       x(:,i) = xin(:,i) + dt*velin(1:ndim,i) + 0.5*dt*dt*forcein(1:ndim,i)           
       vel(:,i) = velin(:,i) + dt*forcein(:,i)
       velin(:,i) = velin(:,i) + 0.5*dt*forcein(:,i)
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
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i) + dt*ddusttogasdtin(i)
          deltav(:,i) = deltavin(:,i) + dt*ddeltavdtin(:,i)
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
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i)
          deltav(:,i)  = deltavin(:,i)
       endif
    else
       vel(:,i) = velin(:,i) + hdt*(force(:,i)) !+forcein(:,i))            
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
       
       if (idust.eq.1) then
          dusttogas(i) = dusttogasin(i) + hdt*(ddusttogasdt(i) + ddusttogasdtin(i))
          deltav(:,i) = deltavin(:,i) + hdt*(ddeltavdt(:,i) + ddeltavdtin(:,i))
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
 if (.not.dtfixed) dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag)
 !if (C_cour*dtav.lt.dt) then
 !   print*,'WARNING: AV controlling timestep: (old)dt = ',dt,' (new)dt = ',C_cour*dtav
 !   dt = C_cour*dtav
 !endif

 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
