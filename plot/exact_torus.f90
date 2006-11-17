! ----------------------------------------------------------------------
!  Plots solution for equilibrium torus of Papaloizou & Pringle
!  strictly valid for the midplane only (uses spherical r)
!
!  Added by D. Price 16.1.06
! ----------------------------------------------------------------------
module torus
  implicit none

contains

subroutine exact_torus(iplot,itorus,Mstar,Rtorus,AA,distortion,gamma,xplot,yplot,ierr)
  implicit none
  integer, intent(in) :: iplot,itorus
  real, intent(in) :: Mstar,Rtorus,AA,gamma,distortion
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  real :: term,densi,rxy
  integer, intent(out) :: ierr
  integer :: i
  real :: ra2
  integer, parameter :: nu = 2
  real, parameter :: atorus = 0.2, currj0 = 1.0
!
! check for errors
!
  ierr = 0
  if (Mstar.le.0.) then
     print*,'error: mass <= 0 in exact_torus'
     ierr = 2
     return
  elseif (Rtorus.lt.0.) then
     print*,'error: rtorus < 0 in exact_torus'
     ierr = 3
     return
  endif
  select case(itorus)

!
!--Tokamak torus (in torus 'r' co-ordinate)
!
  case(2)
  if (nu.le.0 .or. nu.gt.2) then
     print*,'error: solution not found for nu value in tokamak torus'
     ierr = 5
     return
  endif
  
  do i=1,size(xplot)
     ra2 = xplot(i)**2/atorus**2
     if (nu.eq.1) then
        term = currj0**2*atorus**2*(1. - ra2)*(7.*ra2**2 - 23.*ra2 + 13.)/96.
     elseif (nu.eq.2) then
        term = currj0**2*atorus**2*(1. - ra2)*(22.*ra2**4 - 113.*ra2**3 + 237.*ra2**2 - 213.*ra2 + 57.)/720.
     endif

     select case(iplot)
     case(1)
     !--density
        yplot(i) = term
     case(2)
     !--pressure
        yplot(i) = term
     case(3)
     !--thermal energy
        yplot(i) = term
     end select
  enddo    
!
!--Papaloizou & Pringle equilibrium torus
!  
  case default
     if ((gamma-1.).le.1e-4) then
        print*,'error: exact solution not valid for isothermal eos'
        ierr = 4
        return
     endif

     do i=1,size(xplot)
        if (iplot.ne.4) then
           !--plot quantities vs spherical r (assume z = 0)
           term = Mstar/(AA*Rtorus)*(gamma-1.)/gamma* &
                (Rtorus/xplot(i) - 0.5*(Rtorus/xplot(i))**2 - 1./(2.*distortion))
        else
           !--plots with z (assume cyl r = Rtorus)
           rxy = sqrt(Rtorus**2 + xplot(i)**2)
           term = Mstar/(AA*Rtorus)*(gamma-1.)/gamma* &
                (Rtorus/rxy - 0.5 - 1./(2.*distortion))
        endif
        if (term.gt.tiny(term)) then
           densi = term**(1./(gamma-1.))
        else
           densi = 0.
        endif
        select case(iplot)
        case(1)
        !--density
           yplot(i) = densi
        case(2,4)
        !--pressure
           yplot(i) = AA*densi**gamma
        case(3)
        !--thermal energy
           yplot(i) = AA/(gamma-1.)*densi**(gamma-1.)
        end select
     enddo
  end select
  return
end subroutine exact_torus

end module torus
