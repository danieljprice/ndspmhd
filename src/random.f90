!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
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
    iseed = MAX(-iseed,1)  ! do not allow iseed = 0
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

!!------------------------------------------------------------------------!!
!!
!! Long period random number generator (see Press et al, Numerical Recipes)
!! 
!! Period is about 2 x 10**18
!! 
!! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of 
!!  endpoints). Call with iseed < 0 to initialise, thereafter do not
!!  alter iseed between calls.
!!
!!------------------------------------------------------------------------!!

real function ran2(iseed)
 implicit none
 integer, parameter :: im1=2147483563, im2=2147483399, &
   imm1=im1-1, ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, &
   ir2=3791, ntab=32,ndiv=1+imm1/ntab
 real, parameter :: am=1./im1, eps=1.2e-7, rnmx=1.-eps
 integer :: iseed,iseed2,j,k,iv(ntab),iy
 SAVE iv,iy,iseed2

 data iseed2/123456789/, iv/ntab*0/, iy/0/
!
!--initialise random sequence
!
 if (iseed.le.0) then
    iseed = max(-iseed,1) ! iseed not zero
    iseed2 = iseed
    do j=ntab+8,1,-1
       k = iseed/iq1
       iseed = ia1*(iseed-k*iq1) - k*ir1
       if (iseed.lt.0) iseed = iseed + im1
       if (j.le.ntab) iv(j) = iseed
    enddo
    iy = iv(1)
 endif
 k = iseed/iq1
 iseed = ia1*(iseed-k*iq1) - k*ir1
 if (iseed.lt.0) iseed = iseed + im1
 k = iseed2/iq2
 iseed2 = ia2*(iseed2-k*iq2) - k*iq2
 if (iseed2.lt.0) iseed2 = iseed2 + im2
 j = 1 + iy/ndiv
 iy = iv(j) - iseed2
 iv(j) = iseed
 if (iy.lt.1) iy = iy + imm1
 ran2 = min(am*iy,rnmx)
 return
 
end function ran2
 
 
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

!!-------------------------------------------------------------------------
!!
!! Quasi Random sobol sequence from numerical recipes
!!
!!-------------------------------------------------------------------------

SUBROUTINE sobseq(n,x)
 INTEGER n,MAXBIT,MAXDIM
 REAL x(*)
 PARAMETER (MAXBIT=30,MAXDIM=6)
 INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT)
 INTEGER iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
 REAL fac
 SAVE ip,mdeg,ix,iv,in,fac
 EQUIVALENCE (iv,iu)
 DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
 DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
 if (n.lt.0) then
    do k=1,MAXDIM
       do j=1,mdeg(k)
          iu(k,j)=iu(k,j)*2**(MAXBIT-j)
       enddo
       do j=mdeg(k)+1,MAXBIT
          ipp=ip(k)
          i=iu(k,j-mdeg(k))
          i=ieor(i,i/2**mdeg(k))
          do l=mdeg(k)-1,1,-1
            if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
            ipp=ipp/2
          enddo
          iu(k,j)=i
       enddo
    enddo
    fac=1./2.**MAXBIT
    in=0
 else
   im=in
   do j=1,MAXBIT
      if(iand(im,1).eq.0)goto 1
      im=im/2
   enddo
   print*,'MAXBIT too small in sobseq'
1  im=(j-1)*MAXDIM
   do k=1,min(n,MAXDIM)
      ix(k)=ieor(ix(k),iv(im+k))
      x(k)=ix(k)*fac
   enddo
   in=in+1
 endif
 return
END
