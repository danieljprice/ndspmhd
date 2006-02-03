!---------------------------------------------------------------------
! Module containing utilities for the general relativistic code
! These are used by conservative2primitive and primitive2conservative
!---------------------------------------------------------------------
module grutils

contains
!
!--subroutine to calculate the diagonal metric from the coordinates
!  also returns the determinant of the metric sqrtg 
!
  subroutine metric_diag(x,gdiag,sqrtg,ndimx,ndimg,metric)
    implicit none
    integer, intent(in) :: ndimx, ndimg
    real, dimension(ndimx), intent(in) :: x
    real, dimension(ndimg), intent(out) :: gdiag
    real, intent(out) :: sqrtg
    character(len=*), intent(in) :: metric
    real :: rr
    
    select case(metric(1:6))
!
!--cylindrical r,phi,z
!
    case('cylrpz')
       gdiag(:) = 1.
       if (ndimg.ge.2) gdiag(2) = x(1)**2
       sqrtg = abs(x(1))
!
!--cylindrical r,z,phi    
!
    case('cylrzp')
       gdiag(:) = 1.
       if (ndimg.ge.3) gdiag(3) = x(1)**2
       sqrtg = abs(x(1))
!
!--spherical r,phi,theta
!
    case('sphrpt')
        gdiag(:) = 1.
        if (ndimg.ge.2) gdiag(2) = x(1)**2
        if (ndimg.ge.3 .and. ndimx.gt.1) gdiag(3) = gdiag(2)*SIN(x(2))
        sqrtg = x(1)**2
        if (ndimx.gt.1) sqrtg = SIN(x(2))*sqrtg
!
!--spherical with log co-ordinate
!
    case('sphlog')
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
!--subroutine to compute source terms (derivatives of metric) for momentum equation
!  (these are 1/(2*dens)* T^\mu\nu d(g_\mu_\nu)/dx^i
!
  subroutine metric_derivs(xi,veli,Bi,rhoi,pri,B2i,sourceterms,ndim,ndimV,metric)
    implicit none
    integer, intent(in) :: ndim, ndimV
    real, dimension(ndim), intent(in) :: xi
    real, dimension(ndimV), intent(in) :: veli, Bi
    real, intent(in) :: rhoi, pri, B2i
    real, dimension(ndimV), intent(out) :: sourceterms
    character(len=*), intent(in) :: metric
    
    select case(metric(1:6))
    !
    !--cylindrical r,phi,z
    !
    case('cylrpz')
       sourceterms(1) = (pri + 0.5*B2i)/rhoi
       if (ndimV.ge.2) sourceterms(1) = sourceterms(1) + xi(1)*(veli(2)**2 - Bi(2)**2)
    !
    !--cylindrical r,z,phi
    !
    case('cylrzp')
       sourceterms(1) = (pri + 0.5*B2i)/rhoi
       if (ndimV.ge.3) sourceterms(1) = sourceterms(1) + xi(1)*(veli(3)**2 - Bi(3)**2)
    !
    !--spherical
    !
    case('sphrpt')
       sourceterms(1) = (pri + 0.5*B2i)/rhoi
       if (ndimV.ge.2) sourceterms(1) = sourceterms(1) + xi(1)*veli(2)**2
       if (ndimV.ge.3) sourceterms(2) = xi(1)*veli(3)**2
    !
    !--spherical polars with logarithmic radial co-ordinate
    !
    case('sphlog')
       sourceterms(:) = 0.!! WORK THESE OUT!!
    !
    !--note that source terms are zero and not allocated for cartesian
    !
    case default
       sourceterms = 0.
       
    end select
  end subroutine metric_derivs
!
!--function to compute the dot product of two vectors with a diagonal metric
!
  real function DOT_PRODUCT_GR(vec1,vec2,gdiag)
    implicit none
    integer :: i
    real, dimension(:) :: vec1
    real, dimension(size(vec1)) :: vec2,gdiag
    
    DOT_PRODUCT_GR = 0.
    do i=1,size(vec1)
       DOT_PRODUCT_GR = DOT_PRODUCT_GR + gdiag(i)*vec1(i)*vec2(i)
    enddo
    
    return
  end function DOT_PRODUCT_GR
  
end module grutils
