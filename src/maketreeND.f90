!!--------------------------------------------------------------------------
!! MaketreeND - v1.0 1/8/03
!!
!! Constructs the tree for the Barnes/Hut tree code
!! The descriptions are given in Barnes & Hut (1989) ApJS
!!
!! This tree code is designed to work in 1,2 and 3 dimensions
!! Each node in the tree has 2^ndim daughter nodes
!!
!! Written by: Daniel Price, Institute of Astronomy, Cambridge UK
!!--------------------------------------------------------------------------

SUBROUTINE maketree(x,ntot)
! USE debug
! USE loguns
 IMPLICIT NONE		! define local variables
 INTEGER, PARAMETER :: maxlevels = 8
 INTEGER, PARAMETER :: maxdigits = 2**(maxlevels+1)-1, ndim = 3
 INTEGER, INTENT(IN) :: ntot

 INTEGER :: i,j,idim,inode,ilevel
 INTEGER, DIMENSION(ndim) :: ixcoord,idaughter_node
 REAL, DIMENSION(ndim) :: xmin,xmax,dxroot,xcoord
 REAL, DIMENSION(ndim,ntot) :: x
!
!--allow for tracing flow
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine maketree'
!
!--find boundaries of particle distribution
!  (this is the size of the root cell)
! 
 DO j=1,ndim
    xmin(j) = MINVAL(x(j,1:ntot))
    xmax(j) = MAXVAL(x(j,1:ntot))
 ENDDO
 dxroot(:) = xmax(:) - xmin(:)		! size of root cell
!
!--loop over all particles
! 
 over_particles: DO i=1,ntot
!
!--we use a binary representation of the particle's coordinates to determine
!  which node(s) the particle should be in
! 
!--find particle co-ordinates scaled to [0,1) in the root cell
   xcoord(:) = (x(:,i) - xmin(:))/dxroot(:)
!--convert this to an integer value
   ixcoord = NINT(xcoord*maxdigits)
!--use the first two bits of this integer to work out which daughter node
!  the particle should belong to

!
!
!
   ilevel = 0
   nparticles(ilevel) = ntot


   DO WHILE (.NOT.empty)

!
!--shuffle the bits of the coordinates to work out which daughter node the
!  particle should belong to
!
   inode = 0
   DO idim=1,ndim
      CALL MVBITS(ixcoord(idim),maxlevels-ilevel,1,inode,ndim-idim)
   ENDDO
   inode = inode + 1
   
   inode = 0
   idaughter_node = IBITS(ixcoord(:),maxlevels-ilevel,1)
   DO idim=1,ndim
      inode = inode + 2**(ndim-idim)*idaughter_node(idim)
   ENDDO
   inode = inode + 1
   
   PRINT*,' xcoord = ',xcoord,ixcoord,maxdigits
   PRINT*,' bits of xcoord = ',BTEST(ixcoord(1),maxlevels)
   PRINT*,' daughter node = ',idaughter_node,inode
   READ*
!
!--check if this node is empty or not
!
   IF (empty(ilevel,inode)) THEN
!--if empty point to this particle and finish 
   
   ELSE
!--if not empty point to a branch and desend tree
   
   ENDIF
 ENDDO over_particles

 RETURN 
END SUBROUTINE maketree
