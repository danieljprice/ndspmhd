! ----------------------------------------------------------------------
!  Plots various analytic density profiles 
!
!  Currently implemented:
!
!   1) Plummer sphere
!   2) Hernquist sphere
!
!  Added by D. Price 5/9/05
! ----------------------------------------------------------------------
module densityprofiles
  implicit none

contains

subroutine exact_densityprofiles(iplot,iprofile,Msphere,rsoft,xplot,yplot,ierr)
  implicit none
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: iplot,iprofile
  real, intent(in) :: Msphere,rsoft 
  real, intent(in), dimension(:) :: xplot
  real, intent(out), dimension(size(xplot)) :: yplot
  integer, intent(out) :: ierr
  integer :: i
!
! check for errors
!
  ierr = 0
  if (Msphere.le.0.) then
     print*,'error: mass <= 0 in exact_densityprofile'
     ierr = 2
     return
  endif
  if (rsoft.lt.0.) then
     print*,'error: rsoft < 0 in exact_densityprofile'
     ierr = 3
     return
  endif

  select case(iprofile)
  case(1)
!
!--Plummer sphere
!
    select case(iplot)
    case(2)
!--potential
       do i=1,size(xplot)
          yplot(i) = -Msphere/sqrt(rsoft**2 + xplot(i)**2)
       enddo
!--force
    case(3)
       do i=1,size(xplot)
          yplot(i) = Msphere*xplot(i)/((rsoft**2 + xplot(i)**2)**1.5)
       enddo
!--density
    case default
       do i=1,size(xplot)
          yplot(i) = 3.*Msphere*rsoft**2/(4.*pi*(rsoft**2 + xplot(i)**2)**2.5)
       enddo
    end select    

  case(2)
!
!--Hernquist model (use tiny to prevent divergences in cusp)
!
    select case(iplot)
    case(2)
!--potential
       do i=1,size(xplot)
          yplot(i) = -Msphere/sqrt(rsoft**2 + xplot(i)**2)
       enddo
!--force
    case(3)
       do i=1,size(xplot)
          yplot(i) = Msphere*xplot(i)/((rsoft**2 + xplot(i)**2)**1.5)
       enddo 
!--density
    case default
       do i=1,size(xplot)
          yplot(i) = Msphere*rsoft/ &
                     ((2.*pi*xplot(i)*(rsoft + xplot(i))**3))
       enddo
    end select

  case default
     ierr = 1  
  end select
    
  return
end subroutine exact_densityprofiles

end module densityprofiles
