!!-----------------------------------------------------------------
!! This subroutine checks the initial boundary locations and
!! adjusts them if necessary to ensure continuity of the grid
!! setup.
!!
!! THIS SUBROUTINE CURRENTLY ONLY IMPLEMENTED FOR ONE DIMENSION
!!
!!-----------------------------------------------------------------

SUBROUTINE adjust_boundaries(xmin,xmax,x,npart)
 USE dimen_mhd
 USE debug
 USE loguns
!
!--define local variables
!     
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: npart 
 REAL, DIMENSION(npart), INTENT(IN) :: x
 REAL, INTENT(INOUT) :: xmin,xmax
 INTEGER :: ipart,i
 REAL :: difx(5),tol
 LOGICAL :: iuniform
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine adjust_boundaries'
!
!--test to see if grid is uniform near boundaries
!         
 difx(:) = 0.
 tol = 1.e-8	! tolerance should really change with resolution
!
!--left boundary
!	 
 iuniform = .true.
 DO i=1,5
    ipart = i
    difx(i) = x(ipart+1)-x(ipart)
    IF (i.GT.1) THEN
       IF (ABS(difx(i)-difx(i-1)).GE.tol) iuniform = .false.
    ENDIF
 ENDDO
 IF (iuniform) THEN
    xmin = x(1) - 0.5*(x(2)-x(1)) 	
    WRITE(iprint,*) 'uniform grid: xmin adjusted to ',xmin 
 ELSE
    WRITE(iprint,*) 'grid non-uniform: xmin not adjusted'
 ENDIF
!
!--right boundary
!	 
 iuniform = .true.
 DO i=1,5
    ipart = npart-5 + i
    difx(i) = x(ipart)-x(ipart-1)	 
    IF (i.GT.1) THEN
       IF (ABS(difx(i)-difx(i-1)).GE.tol) iuniform = .false.
    ENDIF
 ENDDO
 IF (iuniform) THEN
    xmax = x(npart) + 0.5*(x(npart) - x(npart-1))
    WRITE(iprint,*) 'uniform grid: xmax adjusted to ',xmax
 ELSE
    WRITE(iprint,*) 'grid non-uniform: xmax not adjusted'
 ENDIF

 RETURN
END SUBROUTINE adjust_boundaries
