!
!--calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!
subroutine derivs
 use loguns, only:iprint
 use options, only:ibound,icty
 implicit none
 logical, parameter :: itiming = .false.
 real :: t1,t2,t3,t4,t5
!
!--allow particles to cross boundary (ie. enforce boundary conditions)
!
 if (itiming) call cpu_time(t1)
 if (ANY(ibound.NE.0)) call boundary	! inflow/outflow/periodic boundary conditions
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
 if (icty.LE.0) call iterate_density
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
