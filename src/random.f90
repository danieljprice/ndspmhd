!!------------------------------------------------------------------------!!
!!
!! Random number generator using the minimal standard generator of 
!!  Park & Miller (1988) + shuffling (see Press et al, Numerical Recipes)
!! 
!! Period is about 10**8
!! 
!! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of 
!!  endpoints). Call with iseed < 0 to initialise, thereafter do not
!!  alter iseed between calls.
!!
!!------------------------------------------------------------------------!!

REAL FUNCTION ran1(iseed)
 IMPLICIT NONE
 INTEGER :: iseed
 INTEGER, PARAMETER :: ia = 16807, im=2147483647, iq = 127773, ir = 2836
 INTEGER, PARAMETER :: ntab = 32, ndiv = 1+(im-1)/ntab
 INTEGER, DIMENSION(ntab) :: iv
 INTEGER :: j,k,iy
 REAL, PARAMETER :: am = 1./im, eps = 1.2e-7, floatmax = 1.-eps
 
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
!
!--initialise
! 
 IF (iseed.LE.0 .OR. iy.EQ.0) THEN
    iseed = MAX(-iseed,1) 	! do not allow iseed = 0
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

!!-------------------------------------------------------------------------
!!
!! Function returns a random number drawn from a Rayleigh distribution
!! P(r) = r*e^(-r^2/(2*s^2))/s^2
!! 
!! Useful for drawing amplitudes from a Gaussian distribution,
!! since the modulus is distributed according to a Rayleigh distribution.
!!
!!-------------------------------------------------------------------------
REAL FUNCTION rayleigh_deviate(iseed)
  IMPLICIT NONE
  INTEGER :: iseed
  REAL :: ran1
  
  rayleigh_deviate = SQRT(-LOG(ran1(iseed)))
  
END FUNCTION rayleigh_deviate
