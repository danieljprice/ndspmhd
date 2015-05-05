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

!
!--calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!
subroutine derivs
 use loguns,        only:iprint
 use options,       only:ibound,icty,ihvar,imhd
 use part,          only:hh,x,npart,ntotal,rho,pmass,uu,dustfrac,deltav,pr,spsound
 use rates,         only:ddustevoldt,dendt,drhodt
 use setup_params,  only:hfact
 use cons2prim,     only:conservative2primitive
 use options,       only:idust,iener,use_sqrtdustfrac
 use dustdiffusion, only:dust_diffusion
 use hterms,        only:gradh
 use part,          only:vel
 use rates,         only:force
 implicit none
 logical, parameter :: itiming = .false.
 real :: t1,t2,t3,t4,t5,sum,sum1,sum2,sum3,sum4,sum_dustm,si
 integer :: i,inext
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
!--calculate forces/rates of change using predicted quantities
!
 if (itiming) call cpu_time(t4)

 call get_rates
 
 if (idust.eq.3) then
    !ddustevoldt = 0.
    call dust_diffusion(npart,ntotal,x,pmass,rho,hh,gradh,dustfrac,ddustevoldt,deltav,vel,pr,uu,spsound,dendt)
    sum = 0.
    sum1 = 0.
    sum2 = 0.
    sum3 = 0.
    sum4 = 0.
    sum_dustm = 0.
    do i=1,npart
       sum = sum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - uu(i)*ddustevoldt(i) &
             + (1.-dustfrac(i))*dendt(i))
       sum1 = sum1 + pmass(i)*(dot_product(vel(:,i),force(:,i)))
       sum2 = sum2 - pmass(i)*uu(i)*ddustevoldt(i)
       sum3 = sum3 + pmass(i)*(1. - dustfrac(i))*dendt(i)
       sum4 = sum4 + pmass(i)*0.5*(1. - 2.*dustfrac(i))*ddustevoldt(i)

       if (use_sqrtdustfrac) then
       !
       !--convert from deps/dt to ds/dt
       !
          si = sqrt(dustfrac(i)*rho(i))
          ddustevoldt(i) = 0.5*rho(i)*ddustevoldt(i) + 0.5*si*drhodt(i)/rho(i)
          sum_dustm = sum_dustm + pmass(i)*(2.*si/rho(i)*ddustevoldt(i) - si**2/rho(i)**2*drhodt(i))
       else
          sum_dustm = sum_dustm + pmass(i)*ddustevoldt(i)
       endif
    enddo
    if (abs(sum) > epsilon(sum) .and. iener >= 1) print*,' sum = ',sum,sum1,sum2,sum3,sum4
    if (abs(sum_dustm) > epsilon(sum)) print*,' ERROR in ddustm/dt = ',sum_dustm
 elseif (idust.eq.4 .or. idust.eq.1) then
    sum_dustm = 0.
    if (use_sqrtdustfrac) then
       do i=1,npart
          !
          !--convert from deps/dt to ds/dt
          !
          si = sqrt(dustfrac(i)*rho(i))
          ddustevoldt(i) = 0.5*rho(i)*ddustevoldt(i) + 0.5*si*drhodt(i)/rho(i)
          sum_dustm = sum_dustm + pmass(i)*(2.*si/rho(i)*ddustevoldt(i) - si**2/rho(i)**2*drhodt(i))
       enddo
    else
       do i=1,npart
          sum_dustm = sum_dustm + pmass(i)*ddustevoldt(i)
       enddo
    endif
    if (abs(sum_dustm) > 1.e-14) print*,' ERROR in ddustm/dt = ',sum_dustm
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
