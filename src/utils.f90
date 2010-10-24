!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2010 Daniel Price                                                   !
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

!-------------------------------------------------------------------
! This file contains some simple utilities used in various places
! throughout the code
!------------------------------------------------------------------

!-------------------------------------------------------------------
! simple routine to take min, max and average of a quantity
!-------------------------------------------------------------------
subroutine minmaxave(x,xmin,xmax,xav,npts)
  implicit none
  integer :: i
  integer, intent(in) :: npts
  real, intent(in), dimension(npts) :: x
  real, intent(out) :: xmin,xmax,xav
  
  xav = 0.
  xmin = huge(xmin)
  xmax = -xmin
  do i=1,npts
     xav = xav + x(i)
     xmin = min(xmin,x(i))
     xmax = max(xmax,x(i))
  enddo
  xav = xav/real(npts)

  return
end subroutine minmaxave

subroutine cross_product3D(veca,vecb,vecc)
 implicit none
 real, dimension(3), intent(in) :: veca,vecb
 real, dimension(3), intent(out) :: vecc
 
 vecc(1) = veca(2)*vecb(3) - veca(3)*vecb(2)
 vecc(2) = veca(3)*vecb(1) - veca(1)*vecb(3)
 vecc(3) = veca(1)*vecb(2) - veca(2)*vecb(1)

end subroutine cross_product3D
