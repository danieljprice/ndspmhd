!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!!-----------------------------------------------------------------------
subroutine external_forces(iexternal_force,xpart,fext,ndim)
  use options, only:geom
  use setup_params, only:pi
  implicit none
  integer, intent(in) :: iexternal_force,ndim
  real, dimension(ndim), intent(in) :: xpart
  real, dimension(ndim), intent(out) :: fext
  real, dimension(ndim) :: dr
  real :: rr2,drr,drr2,rcyl2,rcyl,rsph,v2onr
  real, parameter :: Rtorus = 1.0, dfac = 1.1
 real, parameter :: sink = 0.25*pi, Asin = 100., Bsin = 2.0

  fext(:) = 0.
  
  select case(iexternal_force)
  case(1)
!
!--toy star force (linear in co-ordinates)
!
     select case(geom(1:6))
     case ('cylrpz','cylrzp','sphrpt','sphrtp') ! cylindrical, spherical
        fext(1) = -xpart(1)
     case ('cartes')
        fext(:) = - xpart(:)
     case default
        stop 'error: external force not implemented'
     end select

  case(2)
!
!--1/r^2 force from central point mass
!
     dr(:) = 0.
     select case(geom(1:6))
     case('sphrpt','sphrtp') ! spherical
        rr2 = xpart(1)**2
        drr = 1./sqrt(rr2)
        drr2 = drr*drr
        dr(1) = xpart(1)*drr
        fext(:) = -dr(:)*drr2
     case('cylrpz') ! cylindrical r-phi-z
        rr2 = xpart(1)**2
        if (ndim.ge.3) rr2 = rr2 + xpart(3)**2
        drr = 1./sqrt(rr2)
        drr2 = drr*drr 
        dr(1) = xpart(1)*drr
        if (ndim.ge.3) dr(3) = xpart(3)*drr
        fext(:) = -dr(:)*drr2
     case('cylrzp') ! cylindrical r-z-phi
!        rr2 = xpart(1)**2 + xpart(2)**2
!        drr = 1./sqrt(rr2)
!        drr2 = drr*drr
!        dr(1) = xpart(1)*drr
!        dr(2) = xpart(2)*drr
!        fext(:) = -dr(:)*drr2
        fext(:) = -1./xpart(1)**2
     case ('cartes') ! cartesian
        rr2 = DOT_PRODUCT(xpart,xpart) 
        drr = 1./sqrt(rr2)
        drr2 = drr*drr
        dr(:) = xpart(:)*drr
        fext(:) = - dr(:)*drr2
     case default
        stop 'error: external force not implemented'
     end select

  case(3:5)
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
      stop 'external force not implemented'
  case(6)

     fext(:) = 0.
!
!--effective potential for equilibrium torus (Stone, Pringle, Begelman)
!  centripedal force balances pressure gradient and gravity of central point mass 
!
     select case(geom(1:6))
     case('cylrpz')
        rcyl = xpart(1)
        rcyl2 = rcyl*rcyl
        if (ndim.eq.3) then
           rsph = sqrt(rcyl2 + xpart(3)*xpart(3))
        else
           rsph = rcyl
        endif
        v2onr = 1./(Rtorus)*(-Rtorus*rcyl/rsph**3 + Rtorus**2/(rcyl2*rcyl))
        fext(1) = v2onr
        !
        !--for 3D need to add vertical component of gravity
        !
        if (ndim.eq.3) then
           fext(3) = -xpart(3)/rsph**3
        endif     
     case('cartes')
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
     case default
        stop 'external force not implemented'
     end select
  case(7)
!
!--sinusoidal potential as in Dobbs, Bonnell etc.
!
    fext(1) = Asin*sink*SIN(sink*(xpart(1) + Bsin))

  case default
     
     stop 'external force not implemented'
     
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
 use options, only:geom
 use setup_params, only:pi
 implicit none
 integer, intent(in) :: iexternal_force,ndim
 real, dimension(ndim), intent(in) :: xpart
 real, intent(out) :: epot
 real, parameter :: sink = 0.25*pi, Asin = 100., Bsin = 2.0
 
 select case(iexternal_force)
 
 case(1)
!
!--toy star force (x^2 potential)
!
    epot = 0.5*DOT_PRODUCT(xpart,xpart)
 
 case(2)
!
!-- 1/r^2 force(1/r potential)
!
    select case(geom)
    case('sphrtp','sphrpt')
       epot = -1./abs(xpart(1))
    case('cylrpz')
       if (ndim.eq.3) then
          epot = -1./sqrt(xpart(1)**2 + xpart(3)**2)
       else
          epot = -1./abs(xpart(1))
       endif
    case('cylrzp')
       epot = -1./sqrt(DOT_PRODUCT(xpart(1:2),xpart(1:2)))
    case default
       epot = -1./sqrt(DOT_PRODUCT(xpart,xpart))
    end select

 case(3) ! potential from n point masses
    epot = 0.

 case(7)
    epot = Asin*COS(sink*(xpart(1) + Bsin))

 case default
    epot = 0.
 end select

 return
end subroutine external_potentials
