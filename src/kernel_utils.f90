!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
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
!-----------------------------------------------------------------------------
! This module contains utility routines for messing around with the SPH kernel
!-----------------------------------------------------------------------------
module kernel_utils
 implicit none
 integer, parameter  :: dp = 4
 real(dp), parameter :: pi = 3.1415926536

contains
 !-----------------------------------------
 ! routine to normalise the kernel table
 ! appropriately in 1, 2 and 3D
 !-----------------------------------------
 subroutine normalise(ikern,wkern,c,radkern2)
  integer, intent(in)    :: ikern
  real(dp),    intent(in)    :: radkern2
  real(dp),    intent(inout) :: wkern(0:ikern)
  real(dp),    intent(out)   :: c(3)
  real(dp) :: dq2table,q2,q,f(3),qprev,dq
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)

  ! trapezoidal rule to get integral under kernel
  f = 0.
  qprev = 0.
  do i=1,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     dq = q - qprev
     f(1) = f(1) + 0.5*(wkern(i) + wkern(i-1))*dq
     f(2) = f(2) + 0.5*(q*wkern(i) + qprev*wkern(i-1))*dq
     f(3) = f(3) + 0.5*(q2*wkern(i) + qprev**2*wkern(i-1))*dq
     qprev = q
  enddo
  f(1) = 2.*f(1)
  f(2) = 2.*pi*f(2)
  f(3) = 4.*pi*f(3)
  !print*,' integral is ',f(:)
  c(:) = 1./f(:)
  !print*,' normalisation is ',c(:)
 
 end subroutine normalise

 !-----------------------------------------------
 ! differentiate the kernel table to get
 ! tables for the derivative and 2nd deriv
 ! We use 2nd-order accurate finite differencing
 ! on non-uniform meshes
 !-----------------------------------------------
 subroutine diff(ikern,wkern,grkern,grgrkern,radkern2)
  integer, intent(in)    :: ikern
  real(dp), intent(in)    :: radkern2
  real(dp), intent(inout) :: wkern(0:ikern)
  real(dp), intent(out)   :: grkern(0:ikern),grgrkern(0:ikern)
  real(dp) :: dq2table,q2,q,qm1,qm2,qm3,qp1,qp2,qp3
  real(dp) :: dq,dqm1,dqm2,dqm3,dqp1,dqp2,a,b,c,d
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)
  !do i=0,ikern
  !   wkern(i) = i*dq2table
  !enddo
  
  ! finite differencing to get kernel derivatives
  do i=0,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     if (i==0 .or. i==1) then
     ! forward diff
        qp1 = sqrt((i+1)*dq2table)
        qp2 = sqrt((i+2)*dq2table)
        qp3 = sqrt((i+3)*dq2table)
        dq = qp1 - q
        dqp1 = qp2 - qp1
        dqp2 = qp3 - qp2
        a = -(2.*dq + dqp1)/(dq*(dq + dqp1))
        b = (dq + dqp1)/(dq*dqp1)
        c = -dq/(dqp1*(dq + dqp1))      
        grkern(i) = a*wkern(i) + b*wkern(i+1) + c*wkern(i+2)
     ! second deriv, forward diff
        a = (6.d0*dq + 4.d0*dqp1 + 2.d0*dqp2)/(dq*(dq + dqp1)*(dq + dqp1 + dqp2))
        b = -(4.d0*(dq + dqp1) + 2.d0*dqp2)/(dq*dqp1*(dqp1 + dqp2))
        c = (4.d0*dq + 2.d0*(dqp1 + dqp2))/((dqp1 + dq)*dqp1*dqp2)
        d = -(4.d0*dq + 2.d0*dqp1)/((dq + dqp1 + dqp2)*(dqp2 + dqp1)*dqp2)
        grgrkern(i) = a*wkern(i) + b*wkern(i+1) + c*wkern(i+2) + d*wkern(i+3)     
     elseif (i==ikern .or. i==ikern-1) then
     !elseif (i > 3) then
     ! backward diff
        qm1 = sqrt((i-1)*dq2table)
        qm2 = sqrt((i-2)*dq2table)
        qm3 = sqrt((i-3)*dq2table)
        dqm1 = q - qm1
        dqm2 = qm1 - qm2
        dqm3 = qm2 - qm3
        a = dqm1/(dqm2*(dqm1+dqm2))
        b = -(dqm1 + dqm2)/(dqm1*dqm2)
        c = (2.d0*dqm1 + dqm2)/(dqm1*(dqm1 + dqm2))
        grkern(i) = a*wkern(i-2) + b*wkern(i-1) + c*wkern(i)
     ! second deriv, backwards
        a = -(4.d0*dqm1 + 2.d0*dqm2)/(dqm3*(dqm3 + dqm2)*(dqm3 + dqm2 + dqm1))
        b = (4.d0*dqm1 + 2.d0*(dqm2 + dqm3))/(dqm3*dqm2*(dqm2 + dqm1))
        c = -(4.d0*(dqm1 + dqm2) + 2.d0*dqm3)/((dqm2 + dqm3)*dqm2*dqm1)
        d = (6.d0*dqm1 + 4.d0*dqm2 + 2.d0*dqm3)/((dqm1 + dqm2 + dqm3)*(dqm1 + dqm2)*dqm1)
        grgrkern(i) = a*wkern(i-3) + b*wkern(i-2) + c*wkern(i-1) + d*wkern(i)
     else
        qp2 = sqrt((i+2)*dq2table)
        qp1 = sqrt((i+1)*dq2table)
        qm1 = sqrt((i-1)*dq2table)
        dqm1 = q - qm1
        dq   = qp1 - q
        dqp1 = qp2 - qp1
     ! 2nd order unequal grid finite diff
        grkern(i) = -dq/(dqm1*(dq+dqm1))*wkern(i-1) &
                   + (dq-dqm1)/(dq*dqm1)*wkern(i) &
                   + dqm1/(dq*(dq + dqm1))*wkern(i+1)
     ! 2nd deriv - BUT THIS IS ONLY FIRST ORDER ACCURATE
        !a = 2./(dqm1*(dq + dqm1))
        !b = -2./(dqm1*dq)
        !c = 2./(dq*(dq + dqm1))
        !grgrkern(i) = a*wkern(i-1) + b*wkern(i) + c*wkern(i+1)
     ! two nodes forward, one back 2nd derivative estimate
        !a = 2.d0*(2.d0*dq + dqp1)/(dqm1*(dqm1 + dq)*(dqm1 + dq + dqp1))
        !b = -2.d0*(2.d0*dq + dqp1 - dqm1)/(dqm1*dq*(dq + dqp1))
        !c =  2.d0*(dq + dqp1 - dqm1)/((dqm1 + dq)*dq*dqp1)
        !d = -2.d0*(dq - dqm1)/((dqm1 + dq + dqp1)*(dq + dqp1)*dqp1)
        !print*,a,b,c,d
        a = (2.*dqp1 + 4.*dq)/((q-qm1)*(qp1-qm1)*(qp2-qm1))
        b = 2.*((qm1 - q) + 2.*(qp1-q) + qp2-qp1)/((qm1 - q)*(qp1-q)*(qp2 - q))
        c = 2.*(-2.*q + qm1 + qp2)/((q - qp1)*(qm1 - qp1)*(qp2 - qp1))
        d = 2.*(-2.*q + qm1 + qp1)/((q - qp2)*(qm1 - qp2)*(qp1 - qp2))
        !print*,a,b,c,d
        grgrkern(i) = a*wkern(i-1) + b*wkern(i) + c*wkern(i+1) + d*wkern(i+2)
     endif
  enddo
 
 end subroutine diff

 !-----------------------------------------------
 ! differentiate the kernel table to get
 ! tables for the derivative and 2nd deriv
 ! We use 2nd-order accurate finite differencing
 ! on non-uniform meshes
 !-----------------------------------------------
 subroutine differentiate(ikern,wfunc,grfunc,radkern2)
  integer, intent(in)    :: ikern
  real(dp),    intent(in)    :: wfunc(0:ikern),radkern2
  real(dp),    intent(out)   :: grfunc(0:ikern)
  real(dp) :: dq2table,q2,q,qm1,qm2,qp1,qp2
  real(dp) :: dq,dqm1,dqm2,dqp1,a,b,c
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)

  ! finite differencing to get kernel derivatives
  do i=0,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     if (i==0 .or. i==1) then
     ! forward diff
        qp1 = sqrt((i+1)*dq2table)
        qp2 = sqrt((i+2)*dq2table)
        dq = qp1 - q
        dqp1 = qp2 - qp1
        a = -(2.*dq + dqp1)/(dq*(dq + dqp1))
        b = (dq + dqp1)/(dq*dqp1)
        c = -dq/(dqp1*(dq + dqp1))
        grfunc(i) = a*wfunc(i) + b*wfunc(i+1) + c*wfunc(i+2)
     elseif (i==ikern .or. i==ikern-1) then
     ! backward diff
        qm1 = sqrt((i-1)*dq2table)
        qm2 = sqrt((i-2)*dq2table)
        dqm2 = qm1 - qm2
        dqm1 = q - qm1
        a = dqm1/(dqm2*(dqm1+dqm2))
        b = -(dqm1 + dqm2)/(dqm1*dqm2)
        c = (2.*dqm1 + dqm2)/(dqm1*(dqm1 + dqm2))
        grfunc(i) = a*wfunc(i-2) + b*wfunc(i-1) + c*wfunc(i)
     else
        qp1 = sqrt((i+1)*dq2table)
        qm1 = sqrt((i-1)*dq2table)
        dqm1 = q - qm1
        dq   = qp1 - q
     ! 2nd order unequal grid finite diff
        grfunc(i) = -dq/(dqm1*(dq+dqm1))*wfunc(i-1) &
                   + (dq-dqm1)/(dq*dqm1)*wfunc(i) &
                   + dqm1/(dq*(dq + dqm1))*wfunc(i+1)
     endif
  enddo
 
 end subroutine differentiate

 !-----------------------------------------------
 ! differentiate the kernel table to get
 ! tables for the derivative and 2nd deriv
 ! this version does standard finite differencing
 ! on the q2 mesh and then chain rule to get
 ! deriv w.r.t. q instead of q2
 !-----------------------------------------------
 subroutine diffu(ikern,wkern,grkern,grgrkern,radkern2)
  integer, intent(in)    :: ikern
  real(dp),    intent(in)  :: wkern(0:ikern),radkern2
  real(dp),    intent(out) :: grkern(0:ikern),grgrkern(0:ikern)
  real(dp) :: dq2table,q2,q
  real(dp) :: ddq2,ddq22,qp1,qm1
  integer :: i
 
  dq2table = radkern2/real(ikern,kind=dp)
  ddq2     = 1.d0/dq2table
  ddq22    = ddq2*ddq2
  ! finite differencing to get kernel derivatives
  do i=0,ikern
     q2 = i*dq2table
     q  = sqrt(q2)
     if (i==0) then
     ! forward diff
        grkern(i) = 0.5*(-3.*wkern(i) + 4.*wkern(i+1) - wkern(i+2))*ddq2
        grgrkern(i) = (2.*wkern(i) - 5.*wkern(i+1) + 4.*wkern(i+2) - wkern(i+3))*ddq22*2.*q + grkern(i)*2.
        grkern(i) = grkern(i)*2.*q
     elseif (i==ikern) then
     ! backward diff
        grkern(i) = 0.5*(wkern(i-2) - 4.*wkern(i-1) + 3.*wkern(i))*ddq2
        grgrkern(i) = (-wkern(i-3) + 4.*wkern(i-2) - 5.*wkern(i-1) + 2.*wkern(i))*ddq22*2.*q + grkern(i)*2.
        grkern(i) = grkern(i)*2.*q
     else
       qp1 = sqrt((i+1)*dq2table)
       qm1 = sqrt((i-1)*dq2table)
     ! centred diff
        grkern(i) = 0.5*(wkern(i+1) - wkern(i-1))*ddq2
        grgrkern(i) = (wkern(i-1) - 2.*wkern(i) + wkern(i+1))*ddq22*2.*q + grkern(i)*2.
        grkern(i) = grkern(i)*2.*q
        if (i.eq.1 .or. i.eq.2 .or. i.eq.3) print*,i,grkern(i)
     endif
  enddo

 end subroutine diffu

 !-----------------------------------------------
 ! Integrate the kernel table to get
 ! the kernel from the kernel derivative
 ! or kernel derivative from 2nd derivative.
 ! Called twice this does opposite of diff
 !-----------------------------------------------
 subroutine integrate(ikern,df,f,radkern2)
  integer,     intent(in)    :: ikern
  real(dp),    intent(in)    :: df(0:ikern),radkern2
  real(dp),    intent(out)   :: f(0:ikern)
  real(dp) :: dq2table,q,qprev,dq
  integer :: i
  
  dq2table = radkern2/real(ikern,kind=dp)

  ! integrate to get kernel gradient
  f(ikern) = 0.
  qprev = sqrt(ikern*dq2table)
  do i=ikern-1,0,-1
     q  = sqrt(i*dq2table)
     dq = q - qprev
     f(i) = f(i+1) + 0.5*(df(i) + df(i+1))*dq
     qprev = q
  enddo

 end subroutine integrate
end module kernel_utils
