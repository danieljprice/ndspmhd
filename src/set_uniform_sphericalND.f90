!----------------------------------------------------------------
!     Set up a uniform density sphere of particles
!     centred on the origin. This subroutine sets up
!     a uniform cube of particles and trims it to be spherical.
!----------------------------------------------------------------

SUBROUTINE set_uniform_spherical(idist,rmax)
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE part
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: idist
 REAL, INTENT(IN) :: rmax

 INTEGER :: i,j,ierr,ntemp 
 INTEGER, DIMENSION(:), ALLOCATABLE :: partlist
 REAL :: rad
 REAL, DIMENSION(ndim) :: xmin,xmax
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(sphdis)'
!
!--check for errors
!
 IF (rmax.LE.0.) STOP 'Error: set_uniform_spherical: rmax <= 0'
!
!--initially set up a uniform density cartesian grid
!
 xmin(:) = -rmax
 xmax(:) = rmax
 CALL set_uniform_cartesian(idist,psep,xmin,xmax,.false.)
!
!--construct list of particles with r < rmax
!  
 ALLOCATE ( partlist(npart), STAT=ierr )
 IF (ierr.NE.0) STOP 'Error allocating memory in uniform_spherical'

 ntemp = 0		! actual number of particles to use
 DO i=1,npart
    rad = SQRT(DOT_PRODUCT(x(:,i),x(:,i)))
    IF (rad.LE.rmax) THEN
       ntemp = ntemp + 1
       partlist(ntemp) = i
    ENDIF   
 ENDDO 
!
!--reorder particles
!
 npart = ntemp
 ntotal = npart
 x(:,1:npart) = x(:,partlist(1:npart))
!
!--reallocate memory to new size of list
!
 CALL alloc(npart)

 RETURN
END SUBROUTINE
