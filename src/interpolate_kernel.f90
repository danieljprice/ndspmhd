!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and derivative given q^2 = (r_a-r_b)^2/h^2
!!
!! must then divide returned W, grad W by h^ndim, h^ndim+1 respectively
!!----------------------------------------------------------------------

SUBROUTINE interpolate_kernel(q2,W,gradW)
 USE kernel
 IMPLICIT NONE
 INTEGER :: index,index1
 REAL, INTENT(IN) :: q2
 REAL, INTENT(OUT) :: W,gradW
 REAL :: dxx,dwdx,dgrwdx
!
!--find nearest index in kernel table
! 
 index = INT(q2*ddq2table)
 index1 = index + 1
 IF (index.GT.ikern) index = ikern
 IF (index1.GT.ikern) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for W, gradW
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
!
!--interpolate for kernel and derivative
!
 W = (wij(index)+ dwdx*dxx)
 gradW = (grwij(index)+ dgrwdx*dxx)
 
END SUBROUTINE interpolate_kernel
