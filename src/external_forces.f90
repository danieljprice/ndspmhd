!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!! Note that the potential energies for these forces are calculated in
!! evwrite and should correspond to the options used in this subroutine
!!
!!-----------------------------------------------------------------------
subroutine external_forces(ipart,iexternal_force)
  use dimen_mhd
  use part
  use rates
  implicit none
  integer, intent(in) :: ipart,iexternal_force
  real, dimension(ndim) :: dr
  real :: rr,rr2,drr2

  select case(iexternal_force)
  case(1)
!
!--toy star force (linear in co-ordinates)
!
     force(1:ndim,ipart) = force(1:ndim,ipart) - x(1:ndim,ipart)
     
  case(2)
!
!--1/r^2 force from central point mass
!
     rr2 = DOT_PRODUCT(x(:,ipart),x(:,ipart)) 
     rr = SQRT(rr2)
     drr2 = 1./rr2
     dr(:) = x(:,ipart)/rr
     force(1:ndim,ipart) = force(1:ndim,ipart) - dr(1:ndim)*drr2

  case(3)
!
!--compute force on gas particle from n point masses
!
!     do i=1,nptmass
        !
        !--work out distance to point mass
        !        
!        dr(:) = xptmass(:,i) - x(:,ipart)
!        rr2 = DOT_PRODUCT(dr,dr)
!        drr2 = 1./rr2
!        rr = SQRT(rr2)
!        dr(:) = dr(:)/rr
        !
        !--calculate gravitational force
        !  (need to calculate fgrav separately for operator splitting)
        !
!        force(1:ndim,ipart) = force(1:ndim,ipart) - dr(1:ndim)*drr2
!     enddo
  end select
end subroutine external_forces
