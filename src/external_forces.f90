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

module externf
 implicit none
 real, parameter :: eps2_soft = 1.e-8
 real, parameter, public :: Rdisc = 5.0, Mstar = 1.0
 
 public :: external_forces, pequil, external_potentials
 
 private

contains
!!-----------------------------------------------------------------------
!!
!! Computes external (body) forces on a particle given its co-ordinates
!!
!!-----------------------------------------------------------------------
subroutine external_forces(iexternal_force,xpart,fext,ndim,ndimV,vpart,hpart,spsound)
  use options, only:ibound
  use eos,     only:gamma,polyk
  use bound,   only:xmax,xmin
  use part,         only:itypegas,itypedust
  use setup_params, only:xlayer,dwidthlayer,Alayercs,Omega0,Omega2,domegadr,pi
  implicit none
  integer, intent(in) :: iexternal_force,ndim,ndimV
  real, dimension(ndim), intent(in) :: xpart
  real, dimension(ndimV), intent(in) :: vpart
  real, intent(in) :: hpart,spsound
  real, dimension(ndimV), intent(out) :: fext
  real, dimension(ndim) :: dr
  real :: rr,rr2,drr2,rcyl2,rcyl,rsph,v2onr,sink
  real :: q2i,betai,expterm
  real, parameter :: Rtorus = 1.0, dfac = 1.1
  real, parameter :: Asin = 100., Bsin = 2.0
  real, parameter :: smoothl = 0.025
  real :: densmid,denszero,densmedium,yi,dens,ddensdy

  fext(:) = 0.
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
     rr2 = DOT_PRODUCT(xpart,xpart) + eps2_soft
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
     fext(1) = 2.*domegadr*Omega2*xpart(1) + 2.*Omega0*vpart(3)
     fext(2) = 0.
     fext(3) = -2.*Omega0*vpart(1)
!     fext(3) = -Omega2*xpart(3) ! vertical stratification
      !print*,'domegadr = ',2.*domegadr*Omega2,3.*Omega2,xpart(1)
      !print*,'fext = ',fext
      !print*,'vpart = ',vpart      
      !stop 'here'
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
    sink = 0.25*pi
    fext(1) = -Asin*sink*SIN(sink*(xpart(1) + Bsin))

  case(8)
!
!--gravity
!
    if (ndim.ge.2) then
       fext(2) = -0.1
       if (ibound(2).eq.0) then
          q2i = abs(xmax(2)-xpart(2))/hpart
          if (q2i.lt.4.) then
             betai = 0.02*spsound**2/xpart(2)
             fext(2) = fext(2) - betai*gradwkernbound(q2i)
          endif
          q2i = abs(xpart(2)-xmin(2))/hpart
          if (q2i.lt.4.) then
             betai = 0.02*spsound**2/xpart(2)
             fext(2) = fext(2) + betai*gradwkernbound(q2i)
          endif
       endif
    endif
  case(9)
!
!--gravity (2)
!
    if (ndim.ge.2) fext(2) = -10.

  case(10)
!
!--potential for relaxation into Kelvin-Helmholtz initial conditions
!
     denszero = 1.
     densmedium = 2.
     densmid = 0.5*(denszero - densmedium)
     yi = xpart(2) - xmin(2)
     if (yi.gt.0.75) then
        expterm = exp(-(yi-0.75)/smoothl)
        dens    = denszero - densmid*expterm
        ddensdy = densmid/smoothl*expterm
     elseif (yi.gt.0.5) then
        expterm = exp(-(0.75-yi)/smoothl)
        dens    = densmedium + densmid*expterm
        ddensdy = densmid/smoothl*expterm
     elseif (yi.gt.0.25) then
        expterm = exp((-yi + 0.25)/smoothl)
        dens    = densmedium + densmid*expterm
        ddensdy = -densmid/smoothl*expterm
     else
        expterm = exp((yi - 0.25)/smoothl)
        dens    = denszero - densmid*expterm
        ddensdy = -densmid/smoothl*expterm
     endif
     fext(2) = polyk*gamma*dens**(gamma-2.)*ddensdy

  case(12)
!
!--vertical linear force to simulate the vertical motion
!--of particles in a disc
!
     fext(1) = 0.
     fext(2) = - xpart(2) !2D pb
     
  case(13)
!
!--vertical cubic force to benchmark the settling of dust particles
!
     fext(1) = 0.
     fext(2) = - xpart(2)**3 !2D pb 
     
  case(14)
!
!--vertical square root force to benchmark the settling of dust particles
!
     fext(1) = 0.
     if ( xpart(2).ge.0.) then
        fext(2) = - xpart(2)**0.5 !2D pb
     else
        fext(2) = (-xpart(2))**0.5 !2D pb
     endif 

  case(15)
!
!--vertical gravity in disc section (always in LAST coordinate: x in 1D, y in 2D, z in 3D)
!
     fext(:) = 0.
     fext(ndim) = -xpart(ndim)*Mstar/sqrt(Rdisc**2 + xpart(ndim)**2)**3

  case default
     
     fext(:) = 0.
     
  end select
  
  return
end subroutine external_forces

!-----------------------
! kernel boundary force
!-----------------------
real function gradwkernbound(q2)
 implicit none
 real, intent(in) :: q2
 real :: q
 
 q = sqrt(q2)
 if (q.lt.2./3.) then
    gradwkernbound = 2./3.
 elseif (q.lt.1.) then
    gradwkernbound = 2.*q - 1.5*q2
 elseif (q.lt.2.) then
    gradwkernbound = 0.5*(2.-q)**2
 else
    gradwkernbound = 0.
 endif

end function gradwkernbound

real function pequil(iexternal_force,xpart,densi)
 use dimen_mhd, only:ndim
 implicit none
 integer, intent(in) :: iexternal_force
 real, dimension(ndim), intent(in) :: xpart
 real, intent(in) :: densi
 
 select case(iexternal_force)
 case(9)
    pequil = -0.5*densi*(xpart(2)-0.5)
 case(8)
    pequil = -0.1*densi*xpart(2)
 case default
    pequil = 0.
 end select

end function pequil

!!-----------------------------------------------------------------------
!!
!! Calculate potential energy for above forces
!! (called from evwrite in adding up total energy)
!!
!!-----------------------------------------------------------------------

subroutine external_potentials(iexternal_force,xpart,epot,ndim)
 use setup_params, only:pi,domegadr,Omega2
 implicit none
 integer, intent(in) :: iexternal_force,ndim
 real, dimension(ndim), intent(in) :: xpart
 real, intent(out) :: epot
 real, parameter :: Asin = 100., Bsin = 2.0
 real :: sink
 
 select case(iexternal_force)
 case(1) ! toy star force (x^2 potential)
    epot = 0.5*DOT_PRODUCT(xpart,xpart)
 case(2) ! 1/r^2 force(1/r potential)
    epot = -1./SQRT(DOT_PRODUCT(xpart,xpart) + eps2_soft)
 case(3) ! potential from n point masses
    epot = 0.
 case(5)
    epot = domegadr*Omega2*xpart(1)*xpart(1)
 case(7)
    sink = 0.25*pi
    epot = Asin*COS(sink*(xpart(1) + Bsin))
 case default
    epot = 0.
 end select

 return
end subroutine external_potentials

end module externf
