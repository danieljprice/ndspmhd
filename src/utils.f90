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
  xmin = 1.e10
  xmax = -1.e10
  do i=1,npts
     xav = xav + x(i)
     xmin = min(xmin,x(i))
     xmax = max(xmax,x(i))
  enddo
  xav = xav/real(npts)

  return
end subroutine minmaxave

!-------------------------------------------------------------------
! small utility to fit a quadratic function 
! given two points (x0,y0), (x1,y1) and a gradient dydx(x1)
! returns the coefficients of y = ax^2 + bx + c
!-------------------------------------------------------------------
subroutine fit_quadratic(x0,x1,y0,y1,dydx1,aa,bb,cc)
 implicit none
 real, intent(in) :: x0,x1,y0,y1,dydx1
 real, intent(out) :: aa,bb,cc

 aa = ((y1 - y0) - dydx1*(x1-x0))/(x1**2 - x0**2 - 2.*x1*(x1-x0))
 bb = dydx1 - 2.*aa*x1
 cc = y0 - aa*x0**2 - bb*x0

!! print 10,aa,bb,cc
10 format(/,'quadratic: y = ',f9.3,' x^2 + ',f9.3,' x + ',f9.3)

 return
end subroutine fit_quadratic
