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
  subroutine metric_derivs(x,vel,rho,pr,sourceterms,ndim,ndimV,npart,metric)
    implicit none
    integer, intent(in) :: ndim, ndimV, npart
    real, dimension(ndim,npart), intent(in) :: x
    real, dimension(ndimV,npart), intent(in) :: vel
    real, dimension(npart), intent(in) :: rho, pr
    real, dimension(ndimV,npart), intent(out) :: sourceterms
    character(len=*), intent(in) :: metric
    integer :: i
    
    select case(metric(1:6))
    !
    !--cylindrical r,phi,z
    !
    case('cylrpz')
       sourceterms(1,:) = pr(:)/rho(:)
       if (ndimV.ge.2) sourceterms(1,:) = sourceterms(1,:) + x(1,:)*vel(2,:)**2
    !
    !--cylindrical r,z,phi
    !
    case('cylrzp')
       sourceterms(1,:) = pr(:)/rho(:)
       if (ndimV.ge.3) sourceterms(1,:) = sourceterms(1,:) + x(1,:)*vel(3,:)**2
    !
    !--spherical
    !
    case('sphrpt')
       do i=1,size(x(:,i))
          sourceterms(1,i) = pr(i)/rho(i)
          if (ndimV.ge.2) sourceterms(1,i) = sourceterms(1,i) + x(1,i)*vel(2,i)**2
          if (ndimV.ge.3) sourceterms(2,i) = x(1,i)*vel(3,i)**2
       enddo
    !
    !--spherical polars with logarithmic radial co-ordinate
    !
    case('sphlog')
       do i=1,npart
          sourceterms(:,i) = 0.!! WORK THESE OUT!!
       enddo
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
