!--------------------------------------------------------------
! global parameters for NDSPMHD
! this specifies the number of spatial & velocity dimensions
!--------------------------------------------------------------

MODULE dimen_mhd
 IMPLICIT NONE
 INTEGER, PARAMETER :: ndim = 2		! Dimensions
 INTEGER, PARAMETER :: ndimB = 2	! Dimensions of Magnetic field variable
 INTEGER, PARAMETER :: ndimV = ndimB	! Dimensions of Velocity variable
 REAL, PARAMETER :: dndim = 1./REAL(ndim)
END MODULE
