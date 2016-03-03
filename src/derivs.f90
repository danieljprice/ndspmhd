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

!
!--calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!
subroutine derivs
 use loguns,        only:iprint
 use options
 use part,          only:hh,x,npart,ntotal,rho,Bevol,pmass,uu,dustfrac,deltav,pr,spsound,ndust
 use rates,         only:dBevoldt,ddustevoldt,dendt,drhodt
 use setup_params,  only:hfact
 use cons2prim,     only:conservative2primitive
 use resistivity,   only:Bdiffusion
 use timestep,      only:dt
 use dustdiffusion, only:dust_diffusion
 use hterms,        only:gradh
 use part,          only:vel,P_Q,rhodust,rhogas
 use rates,         only:force,gradpsi
 use get_quantum,   only:get_quantum_pressure
 use aniso_diffusion, only:get_aniso_diffusion
 use timestep,        only:dtcourant,dtforce,dtdrag,dtvisc
 use getgrad,         only:get_gradient
 use getdiv,          only:get_divergence
 use utils,           only:delta_fn
 use bound,           only:ireal
 implicit none
 logical, parameter :: itiming = .false.
 integer :: i,j,inext
 real :: t1,t2,t3,t4,t5,sum,sum1,sum2,sum3,sum4,sum_dustm,dustm,si(ndust)
 real :: k_iso, k_par
 real, dimension(3) :: Bvec
 real, dimension(ndim) :: tmp
 real, dimension(ndim,ndim) :: k_tensor
 real, save :: dustm_prev = 0.
!
!--allow particles to cross boundary (ie. enforce boundary conditions)
!
 if (itiming) call cpu_time(t1)
 
 if (ihvar.EQ.5) then
    inext = int(2.*hfact)
    !!print*,'inext = ',inext
    do i=inext+1,npart-inext
       hh(i) = 0.5*max(abs(x(1,i+inext)-x(1,i))+abs(x(1,i+inext-1)-x(1,i)), &
                       abs(x(1,i-inext)-x(1,i))+abs(x(1,i-inext+1)-x(1,i)))
    enddo
    do i=1,inext
       hh(i) = hh(inext)
    enddo
    do i=npart-inext+1,npart
       hh(i) = hh(npart-inext+1)
    enddo
 endif
 if (ANY(ibound.NE.0)) call boundary ! inflow/outflow/periodic boundary conditions
!
!--set ghost particles if ghost boundaries are used
!         
 if (ANY(ibound.GE.2)) call set_ghost_particles
!
!--call link list to find neighbours
!
 call set_linklist
!
!--calculate density by direct summation
!
 if (itiming) call cpu_time(t2)
 if (imhd.eq.5) then
    imhd = 55
    call iterate_density
    imhd = 5
 else
    if (icty.LE.0) call iterate_density
 endif
!
!--calculate primitive variables from conservative variables
!
 if (itiming) call cpu_time(t3)
 call conservative2primitive
!
!--calculate quantum pressure term (if required)
!
 if (iquantum /= 0 .and. iquantum.ne.2 .and. iquantum.ne.3) then
    call get_quantum_pressure(iquantum,npart,x,pmass,rho,hh,P_Q)
 endif
!
!--calculate forces/rates of change using predicted quantities
!
 if (itiming) call cpu_time(t4)

 if (idiffuse > 0) then
    !
    ! use this for performing tests of anisotropic diffusion
    !
    if (iener /= 2) stop 'anisotropic diffusion assumes iener=2'
    Bvec = (/0.,0.,1./)
    k_iso = 0.0!2
    k_par = 0.01
    do i=1,ndim
       do j=1,ndim
          k_tensor(i,j) = k_iso*delta_fn(i,j) + k_par*Bvec(i)*Bvec(j)
       enddo
    enddo
    force(:,:) = 0.

    if (idiffuse==2) then
       ! take gradient with differencing operator
       call get_gradient(1,npart,x,pmass,rho,hh,uu,gradpsi)

       ! copy to ghosts
       do i=npart+1,ntotal
          j = ireal(i)
          if (j > 0) gradpsi(:,i) = gradpsi(:,j)
       enddo

       ! dot with diffusion tensor
       do i=1,ntotal
          do j=1,ndim
             tmp(j) = -dot_product(k_tensor(j,:),gradpsi(:,i))
          enddo
          gradpsi(:,i) = -tmp(:)
       enddo

       ! take divergence with symmetric operator
       call get_divergence(2,npart,x,pmass,rho,hh,gradpsi,dendt)
    elseif (idiffuse==1) then
       call get_aniso_diffusion(npart,x,pmass,rho,hh,uu,dendt,k_tensor)
    endif
    dtcourant = huge(dtcourant)
    dtforce   = huge(dtforce)
    dtdrag    = huge(dtdrag)
    dtvisc    = huge(dtvisc)
    do i=1,npart
       dtcourant = min(dtcourant,hh(i)**2/(k_iso + k_par),0.1)
    enddo
 else
    call get_rates

    if (imhd.eq.11 .and. iresist.eq.2 .and. etamhd.gt.0.) then
       call Bdiffusion(npart,x,pmass,rho,hh,Bevol,dBevoldt,dt)
    endif
 endif
 
 if (idust.eq.3) then
    !ddustevoldt = 0.
    call dust_diffusion(npart,ntotal,x,pmass,rho,hh,gradh,dustfrac,ddustevoldt,deltav,vel,pr,uu,spsound,dendt)
    sum = 0.
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum_dustm = 0.
    dustm = 0.
    do i=1,npart
       sum = sum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - uu(i)*ddustevoldt(1,i) &
             + (1.-dustfrac(1,i))*dendt(i))
       sum1 = sum1 + pmass(i)*(dot_product(vel(:,i),force(:,i)))
       sum2 = sum2 - pmass(i)*uu(i)*ddustevoldt(1,i)
       sum3 = sum3 + pmass(i)*(1. - dustfrac(1,i))*dendt(i)
       sum4 = sum4 + pmass(i)*0.5*(1. - 2.*dustfrac(1,i))*ddustevoldt(1,i)

       select case(idustevol)
       case(3)
          ddustevoldt(:,i) = 2.*ddustevoldt(:,i)
          sum_dustm = sum_dustm + pmass(i)*0.5*ddustevoldt(1,i)
       case(2)
          ddustevoldt(:,i) = 0.5*rho(i)*ddustevoldt(:,i)
          if (use_smoothed_rhodust) then
             si(:) = sqrt(rhodust(:,i)*rhogas(i))
          else
             si(:) = sqrt(rho(i)**2*(dustfrac(:,i)*(1. - dustfrac(:,i))))
          endif
          sum_dustm = sum_dustm + pmass(i)*si(1)/rho(i)*ddustevoldt(1,i)
       case(1)
       !
       !--convert from deps/dt to ds/dt
       !
          si(:) = sqrt(dustfrac(:,i)*rho(i))
          ddustevoldt(:,i) = 0.5*rho(i)*ddustevoldt(:,i) + 0.5*si(:)*drhodt(i)/rho(i)
          sum_dustm = sum_dustm + pmass(i)*(2.*si(1)/rho(i)*ddustevoldt(1,i) - si(1)**2/rho(i)**2*drhodt(i))
       case default
          sum_dustm = sum_dustm + pmass(i)*ddustevoldt(1,i)
       end select
       dustm = dustm + pmass(i)*dustfrac(1,i)
    enddo
    if (abs(sum) > epsilon(sum) .and. iener >= 1) print*,' sum = ',sum,sum1,sum2,sum3,sum4
    if (abs(sum_dustm) > epsilon(sum)) print*,' ERROR in ddustm/dt = ',sum_dustm
    print*,' dustm = ',dustm
 elseif (idust.eq.4 .or. idust.eq.1) then
    sum_dustm = 0.
    dustm = 0.
    select case(idustevol)
    case(3)
       do i=1,npart
          ddustevoldt(:,i) = 2.*ddustevoldt(:,i)
          sum_dustm = sum_dustm + pmass(i)*0.5*ddustevoldt(1,i)
       enddo
    case(2)
       do i=1,npart
          ddustevoldt(:,i) = 0.5*ddustevoldt(:,i)
          if (use_smoothed_rhodust) then
             si(:) = sqrt(rhodust(:,i)*rhogas(i))
          else
             si(1) = rho(i)*sqrt((dustfrac(1,i)*(1. - dustfrac(1,i))))
          endif
          sum_dustm = sum_dustm + 2.*pmass(i)*si(1)/rho(i)*ddustevoldt(1,i)
          dustm = dustm + pmass(i)*dustfrac(1,i)
       enddo
    case(1)
       do i=1,npart
          !
          !--convert from deps/dt to ds/dt
          !
          if (use_smoothed_rhodust) then
             si(:) = sqrt(rhodust(:,i))
          else
             si(:) = sqrt(rho(i)*dustfrac(:,i))
          endif
          ddustevoldt(:,i) = 0.5*rho(i)*ddustevoldt(:,i) + 0.5*si(:)*drhodt(i)/rho(i)
          sum_dustm = sum_dustm + pmass(i)*(2.*si(1)/rho(i)*ddustevoldt(1,i) - si(1)**2/rho(i)**2*drhodt(i))
          dustm = dustm + pmass(i)*dustfrac(1,i)
       enddo
    case default
       do i=1,npart
          sum_dustm = sum_dustm + pmass(i)*ddustevoldt(1,i)
          dustm = dustm + pmass(i)*dustfrac(1,i)
       enddo
    end select
    print*,' dustm = ',dustm,(dustm-dustm_prev)*pmass(1),sum_dustm
    dustm_prev = dustm
    if (abs(sum_dustm) > epsilon(sum)) print*,' ERROR in ddustm/dt = ',sum_dustm
 endif

 if (itiming) then
    call cpu_time(t5)
    write(iprint,"(50('-'))") 
    write(iprint,*) 'time for dens    = ',t3-t2,'s'
    write(iprint,*) 'time for rates   = ',t5-t4,'s'
    write(iprint,*) '        others   = ',(t2-t1)+(t4-t3),'s'
    write(iprint,*) 'total for derivs = ',t5-t1,'s'
    write(iprint,"(50('-'))") 
 endif
 
 return
end subroutine derivs
