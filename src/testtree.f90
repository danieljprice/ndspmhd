!!
!! Testtree - tests the maketree subroutine
!! sets up a random distribution of particles and calls tree construction
!!
PROGRAM testtree
 IMPLICIT NONE
 INTEGER, PARAMETER :: ndim = 3
 INTEGER, PARAMETER :: npart = 500
 INTEGER :: i,j,iseed
 REAL, DIMENSION(ndim,npart) :: x
 REAL :: ran1
!
!--setup particles
! 
 iseed = -2356
 x(1,1) = ran1(iseed)
 DO i=1,npart
    DO j = 1,ndim
       x(j,i) = ran1(iseed)
    ENDDO   
 ENDDO
!
!--construct the tree code
!
 CALL maketree(x,npart)

END PROGRAM testtree

!!------------------------------------------------------------------------!!
!!
!! Random number generator from numerical recipes
!! 
!! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of 
!!  endpoints). Call with iseed < 0 to initialise, thereafter do not
!!  alter iseed between calls.
!!
!!------------------------------------------------------------------------!!

FUNCTION ran1(iseed)
 IMPLICIT NONE
 INTEGER :: iseed
 INTEGER, PARAMETER :: ia = 16807, im=2147483647, iq = 127773, ir = 2836
 INTEGER, PARAMETER :: ntab = 32, ndiv = 1+(im-1)/ntab
 INTEGER, DIMENSION(ntab) :: iv
 INTEGER :: j,k,iy
 REAL :: ran1
 REAL, PARAMETER :: am = 1./im, eps = 1.2e-7, floatmax = 1.-eps
 
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
!
!--initialise
! 
 IF (iseed.LE.0 .OR. iy.EQ.0) THEN
    iseed = MAX(-iseed,1)       ! do not allow iseed = 0
    DO j = ntab+8,1,-1
       k = iseed/iq
       iseed = ia*(iseed-k*iq) - ir*k
       IF (iseed.LT.0) iseed = iseed + im
       IF (j.LE.ntab) iv(j) = iseed
    ENDDO
    iy = iv(1)
 ENDIF
!
!--generate random number
!
 k = iseed/iq
 iseed = ia*(iseed-k*iq) - ir*k
 IF (iseed.LT.0) iseed = iseed + im
 j = 1 + iy/ndiv
 iy = iv(j)
 iv(j) = iseed
 ran1 = MIN(AM*iy,floatmax)
 
 RETURN
END FUNCTION ran1
