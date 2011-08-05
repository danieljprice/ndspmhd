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

subroutine curl3D_epsijk(gradAvec,curlA)
 implicit none
 real, dimension(3,3), intent(in) :: gradAvec
 real, dimension(3), intent(out) :: curlA
 
 curlA(1) = gradAvec(2,3) - gradAvec(3,2)
 curlA(2) = gradAvec(3,1) - gradAvec(1,3)
 curlA(3) = gradAvec(1,2) - gradAvec(2,1)

end subroutine curl3D_epsijk
