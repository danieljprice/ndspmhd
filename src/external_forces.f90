!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!!-----------------------------------------------------------------------
subroutine external_forces(iexternal_force,xpart,fext,ndim,ndimV,vpart)
  use setup_params, only:xlayer,dwidthlayer,Alayercs,Omega0,Omega2,domegadr,pi
  implicit none
  integer, intent(in) :: iexternal_force,ndim
  real, dimension(ndim), intent(in) :: xpart
  real, dimension(ndimV), intent(in) :: vpart
  real, dimension(ndim), intent(out) :: fext
  real, dimension(ndim) :: dr
  real, intent(in) :: vphi
  real :: rr,rr2,drr2,rcyl2,rcyl,rsph,v2onr,drcyl(2)
  real, parameter :: Rtorus = 1.0, dfac = 1.1
  real, parameter :: sink = 0.25*pi, Asin = 100., Bsin = 2.0

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
!--this is for the 2D cartesian shearing box
!
     fext(1) = 2.*domegadr*Omega2*xpart(1) + 2.*Omega0*vpart(2)
     fext(2) = -2.*Omega0*vpart(1)
!     fext(3) = -Omega2*xpart(3) ! vertical stratification
  case(6)
     fext(:) = 0.
!
!--effective potential for equilibrium torus (Stone, Pringle, Begelman)
!  centripedal force balances pressure gradient and gravity of central point mass 
!
     rcyl2 = DOT_PRODUCT(xpart(1:2),xpart(1:2))
     rcyl = SQRT(rcyl2)
     if (ndim.eq.3) then
        rsph = sqrt(rcyl2 + xpart(3)*xpart(3))
     else
        rsph = rcyl
     endif
     v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl))
     fext(1:2) = v2onr*xpart(1:2)/rcyl
!
!--for 3D need to add vertical component of gravity
!
     if (ndim.eq.3) then
        fext(3) = -xpart(3)/rsph**3
     endif
     
  case(7)
!
!--sinusoidal potential as in Dobbs, Bonnell etc.
!
    fext(1) = -Asin*sink*SIN(sink*(xpart(1) + Bsin))

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

subroutine external_potentials(iexternal_force,xpart,epot,ndim,ndimV,vpart)
 implicit none
 integer, intent(in) :: iexternal_force,ndim,ndimV
 real, dimension(ndim), intent(in) :: xpart
 real, intent(out) :: epot
 real, dimension(ndimV), intent(in) :: vpartc
 real, parameter :: sink = 0.25*pi, Asin = 100., Bsin = 2.0
 
 select case(iexternal_force)
 case(1) ! toy star force (x^2 potential)
    epot = 0.5*DOT_PRODUCT(xpart,xpart)
 case(2) ! 1/r^2 force(1/r potential)
    epot = -1./SQRT(DOT_PRODUCT(xpart,xpart))
 case(3) ! potential from n point masses
    epot = 0.
 case(5)
    epot = domegadr*Omega2*xpart(1)*xpart(1)
 case(7)
    epot = Asin*COS(sink*(xpart(1) + Bsin))
 case default
    epot = 0.
 end select

 return
end subroutine external_potentials
