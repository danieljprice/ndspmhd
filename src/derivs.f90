!
!--calculates derivatives of all particle quantities
!  (wrapper for call to density and rates, calls neighbours etc first)
!
subroutine derivs
 use options, only:ibound,icty
 implicit none
!
!--allow particles to cross boundary (ie. enforce boundary conditions)
!
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
 if (icty.LE.0) call iterate_density
!
!--calculate primitive variables from conservative variables
!   
 call conservative2primitive
!
!--calculate forces/rates of change using predicted quantities
!
 call get_rates

 return
end subroutine derivs
