!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!!-----------------------------------------------------------------------
subroutine external_forces(iexternal_force,xpart,fext,ndim)
  implicit none
  integer, intent(in) :: iexternal_force,ndim
  real, dimension(ndim), intent(in) :: xpart
  real, dimension(ndim), intent(out) :: fext
  real, dimension(ndim) :: dr
  real :: rr,rr2,drr2

  select case(iexternal_force)
  case(1)
!
!--toy star force (linear in co-ordinates)
!
     fext(:) = - xpart(:)
     
  case(2)
!
!--1/r^2 force from central point mass
!
     rr2 = DOT_PRODUCT(xpart,xpart) 
     rr = SQRT(rr2)
     drr2 = 1./rr2
     dr(:) = xpart(:)/rr
     fext(:) = - dr(:)*drr2

  case(3)
!
!--compute force on gas particle from n point masses
!
     fext = 0.
!     do i=1,nptmass
        !
        !--work out distance to point mass
        !        
!        dr(:) = xptmass(:,i) - xpart(:)
!        rr2 = DOT_PRODUCT(dr,dr)
!        drr2 = 1./rr2
!        rr = SQRT(rr2)
!        dr(:) = dr(:)/rr
        !
        !--calculate gravitational force
        !
!        fext(:) = fext(:) - dr(:)*drr2
!     enddo
  case default
     
     fext(:) = 0.
     
  end select
  
  return
end subroutine external_forces

!!-----------------------------------------------------------------------
!!
!! Calculate potential energy for above forces
!! (called from evwrite in adding up total energy)
!!
!!-----------------------------------------------------------------------

subroutine external_potentials(iexternal_force,xpart,epot,ndim)
 implicit none
 integer, intent(in) :: iexternal_force,ndim
 real, dimension(ndim), intent(in) :: xpart
 real, intent(out) :: epot
 
 select case(iexternal_force)
 case(1) ! toy star force (x^2 potential)
    epot = 0.5*DOT_PRODUCT(xpart,xpart)
 case(2) ! 1/r^2 force(1/r potential)
    epot = 1./SQRT(DOT_PRODUCT(xpart,xpart))
 case(3) ! potential from n point masses
    epot = 0.
 case default
    epot = 0.
 end select

 return
end subroutine external_potentials
