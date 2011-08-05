!
!--calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!
subroutine derivs
 use loguns, only:iprint
 use options, only:ibound,icty,ihvar,imhd,isplitpart
 use part, only:hh,x,npart,rho,Bevol,pmass
 use rates, only:dBevoldt
 use setup_params, only:hfact
 use cons2prim, only:conservative2primitive
 use resistivity, only:Bdiffusion
 use timestep, only:dt,C_Cour,C_force,dtcourant,dtforce
 use options, only:iresist,etamhd
 implicit none
 logical, parameter :: itiming = .false.
 real :: t1,t2,t3,t4,t5
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
 
 if (imhd.eq.11 .and. iresist.eq.2 .and. etamhd.gt.0.) then
    !print*,' in derivs, dt = ',dt,min(C_force*dtforce,C_cour*dtcourant)
    call Bdiffusion(npart,x,pmass,rho,hh,Bevol,dBevoldt,min(C_force*dtforce,C_cour*dtcourant))
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
