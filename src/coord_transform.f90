!-----------------------------------------------------------------
! Subroutine to transform between different co-ordinate systems
! (e.g. from cartesian to cylindrical polar and vice versa)
!
! xin(ndimin)   : input co-ordinates, in ndimin dimensions
! itypein       : input co-ordinate type
!
! xout(ndimout) : output co-ordinates, in ndimout dimensions
! itypeout      : output co-ordinate type
!
! itype must be one of the following:
!  itype = 1    : cartesian (default)
!  itype = 2    : cylindrical
!  itype = 3    : spherical
!
! Currently handles:
!
!  cartesian -> cylindrical, spherical polar
!  cylindrical -> cartesian
!  spherical polar -> cartesian
!
!-----------------------------------------------------------------
subroutine coord_transform(xin,ndimin,itypein,xout,ndimout,itypeout)
  implicit none
  integer, intent(in) :: ndimin,ndimout,itypein,itypeout
  real, intent(in), dimension(ndimin) :: xin
  real, intent(out), dimension(ndimout) :: xout
!
!--check for errors in input
!
  if (itypeout.eq.itypein) then
     xout(1:ndimout) = xin(1:ndimout)
     return
  elseif (ndimin.lt.1.or.ndimin.gt.3) then
     print*,'Error: coord transform: invalid number of dimensions on input'
     return
  elseif (ndimout.lt.1.or.ndimout.gt.3) then
     print*,'Error: coord transform: invalid number of dimensions on output'
     return
  elseif (ndimout.gt.ndimin) then
     print*,'Error: coord transform: ndimout must be <= ndimin'
     return
  elseif (abs(xin(1)).lt.1e-8 .and. ndimout.ge.2 .and. &
       (itypeout.eq.2 .or. itypeout.eq.3)) then
     print*,'Warning: coord transform: r=0 on input: cannot return angle'
     xout(1:ndimout) = xin(1:ndimout)
     return
  endif
!
!--now do transformation
!
  select case(itypein)
!
!--input is cylindrical polars
!
  case(2)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        if (itypeout.ne.1) print*,'warning: using default cartesian output'
        if (ndimout.eq.1) then
           xout(1) = xin(1)
        else  ! r,phi,z -> x,y,z
           xout(1) = xin(1)*COS(xin(2))
           xout(2) = xin(1)*SIN(xin(2))
           if (ndimout.gt.2) xout(3) = xin(3)
        endif
     end select
!
!--input is spherical polars
!
  case(3)
     select case(itypeout)
     case default
        !        
        ! output is cartesian (default)
        !
        if (itypeout.ne.1) print*,'warning: using default cartesian output'
        select case(ndimout)
           case(1) ! r -> x
              xout(1) = xin(1)
           case(2) ! r,phi -> x,y
              xout(1) = xin(1)*COS(xin(2))
              xout(2) = xin(1)*SIN(xin(2))
           case(3) ! r,theta,phi -> x,y,z
              xout(1) = xin(1)*SIN(xin(2))*COS(xin(3))
              xout(2) = xin(1)*SIN(xin(2))*SIN(xin(3))
              xout(3) = xin(1)*COS(xin(2))
        end select
     end select
!
!--input is cartesian co-ordinates
!
  case default
     select case(itypeout)
     case(2)
        !
        !--output is cylindrical
        !
        if (ndimin.eq.1) then
           xout(1) = abs(xin(1))   ! cylindrical r
        else
           xout(1) = SQRT(DOT_PRODUCT(xin(1:2),xin(1:2)))
           if (ndimout.ge.2) xout(2) = ATAN(xin(2)/xin(1)) ! phi
           if (ndimout.eq.3) xout(3) = xin(3)              ! z
        endif
     case(3)
        !
        ! output is spherical
        !
        xout(1) = SQRT(DOT_PRODUCT(xin,xin))! r  
        if (ndimout.eq.2) then              ! 2D spherical returns r,phi
           xout(2) = ATAN(xin(2)/xin(1))    ! phi = ATAN(y/x)
        elseif (ndimout.eq.3) then
           xout(2) = ACOS(xin(3)/xout(1))   ! theta = ACOS(z/r)
           xout(3) = ATAN(xin(2)/xin(1))    ! phi = ATAN(y/x)
        endif
     case default
        !
        ! just copy
        !
	xout(1:ndimout) = xin(1:ndimout)
     end select
  end select

  return
end subroutine coord_transform
