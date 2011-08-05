module errors
!
! Fortran 90 module containing subroutine for calculation of various error norms
! input is tabulated exact solution (xexact, yexact) and tabulated approximate
! solution (xpts, ypts). Returns residuals at each of the xpts 
! using linear interpolation from and the xexact and yexact points and
! L1, L2 and L_inf error norms for the points supplied. Only uses points that
! are within range of xexact.
!
! (c) 2006 Daniel Price, University of Exeter dprice@astro.ex.ac.uk
!
!
  implicit none
  public :: calculate_errors

  private
 
contains

  subroutine calculate_errors(xexact,yexact,xpts,ypts,residual,errL1,errL2,errLinf)
    implicit none
    real, dimension(:), intent(in) :: xexact,yexact,xpts,ypts
    real, dimension(size(xpts)), intent(out) :: residual
    real, intent(out) :: errL1,errL2,errLinf
    integer :: i,j,npart,iused
    real :: xi,dy,dx,yexacti,err1,ymax

    errL1 = 0.
    errL2 = 0.
    errLinf = 0.
    residual = 0.
    npart = size(xpts)
    iused = 0
    ymax = -huge(ymax)

    do i=1,npart
       xi = xpts(i)
       yexacti = 0.
       !
       !--find nearest point in exact solution table
       !
       do j=1,size(xexact)-1
          if (xexact(j).le.xi .and. xexact(j+1).gt.xi) then
             if (abs(residual(i)).gt.tiny(residual)) print*,'already used ',i,residual(i)
             !--linear interpolation from tabulated exact solution
             dy = yexact(j+1) - yexact(j)
             dx = xexact(j+1) - xexact(j)
             if (dx.gt.0.) then
                yexacti = yexact(j) + dy/dx*(xi - xexact(j))
                residual(i) = ypts(i) - yexacti
             elseif (dy.gt.0.) then
                yexacti = yexact(j)
                residual(i) = ypts(i) - yexacti
             else
                print "(a)",'error in residual calculation'
                residual(i) = 0.
             endif
             iused = iused + 1
             ymax = max(ymax,abs(yexacti))
          endif
       enddo
       err1 = abs(residual(i))
       errL1 = errL1 + err1
       errL2 = errL2 + err1**2
       errLinf = max(errLinf,err1)
       if (yexacti.gt.tiny(yexacti)) residual(i) = residual(i)/abs(yexacti)
    enddo
    !
    !--normalise errors (use maximum y value)
    !
    if (ymax.gt.tiny(ymax)) then
       errL1 = errL1/(npart*ymax)
       errL2 = sqrt(errL2/(npart*ymax**2))
       errLinf = errLinf/ymax
    else
       print "(a)",'error normalising errors'
       errL1 = 0.
       errL2 = 0.
       errLinf = 0.
    endif

    if (iused.ne.npart) print*,'errors calculated using ',iused,' of ',npart, 'particles'

    return
  end subroutine calculate_errors

end module errors
