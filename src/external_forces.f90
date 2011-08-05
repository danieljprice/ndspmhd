!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!!-----------------------------------------------------------------------
subroutine external_forces(iexternal_force,xpart,fext,ndim,vphi)
  use setup_params, only:xlayer,dwidthlayer,Alayercs,Omega,Omega2
  implicit none
  integer, intent(in) :: iexternal_force,ndim
  real, dimension(ndim), intent(in) :: xpart
  real, dimension(ndim), intent(out) :: fext
  real, dimension(ndim) :: dr
  real, intent(in) :: vphi
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
  case(4)
!
!--boundary layer for acoustic-shock problem
!
     fext(:) = 0.
     fext(1) = Alayercs*dwidthlayer*1./(COSH((xpart(1)-xlayer)*dwidthlayer))**2

  case(5)
!
!--external forces for the 2D MRI problem (ie. central point mass, centrifugal forces)
!
!
!--1/r^2 force from central point mass (r is relative to grid centre)
!
!     rr = Rcentre + xpart(1)
!     fext(1) = -1./rr**2
!
!--add centrifugal force (assume keplerian rotation)
!
!     fext(1) = 1./Rcentre**2
!
!--coriolis force (vphi should be set such that the forces balance initially)
!
     fext(1) = 3.*Omega2*xpart(1) + 2.*Omega*vphi
     fext(:) = 0.
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
