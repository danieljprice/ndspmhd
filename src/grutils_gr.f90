!---------------------------------------------------------------------
! Module containing utilities for the general relativistic code
! These are used by conservative2primitive and primitive2conservative
!---------------------------------------------------------------------
module grutils

contains
!
!--subroutine to calculate the diagonal metric from the coordinates
!  also returns the determinant of the metric sqrtg 
!  (note that derivatives still need to be known in conservative2primitive)
!
  subroutine metric_diag(x,gdiag,sqrtg,ndimx,ndimg,imetric)
    implicit none
    integer :: ndimx, ndimg, imetric
    real, dimension(ndimx) :: x
    real, dimension(ndimg) :: gdiag
    real :: sqrtg, rr
    
    select case(imetric)
!
!--cylindrical
!
    case(1)
       gdiag(:) = 1.
       if (ndimg.ge.2) gdiag(2) = x(1)**2
       sqrtg = abs(x(1))
!
!--spherical
!
    case(2)
        gdiag(:) = 1.
        if (ndimg.ge.2) gdiag(2) = x(1)**2
        if (ndimg.ge.3 .and. ndimx.gt.1) gdiag(3) = gdiag(2)*SIN(x(2))
        sqrtg = x(1)**2
        if (ndimx.gt.1) sqrtg = SIN(x(2))*sqrtg
!
!--spherical with log co-ordinate
!
    case(3)
       gdiag(:) = 1.
       rr = exp(x(1))
       if (ndimx.eq.1) then
          sqrtg = rr**3
       else
          sqrtg = SIN(x(2))*rr**3
       endif
!
!--cartesian
!
    case default
       gdiag = 1.
       sqrtg = 1.
    end select
    
  end subroutine metric_diag

!
!--function to compute the dot product of two vectors with a diagonal metric
!
  real function DOT_PRODUCT_GR(vec1,vec2,gdiag,ndim)
    implicit none
    integer :: i,ndim
    real, dimension(ndim) :: vec1,vec2,gdiag
    
    DOT_PRODUCT_GR = 0.
    do i=1,ndim
       DOT_PRODUCT_GR = DOT_PRODUCT_GR + gdiag(i)*vec1(i)*vec2(i)
    enddo
    
    return
  end function DOT_PRODUCT_GR
  
end module grutils
