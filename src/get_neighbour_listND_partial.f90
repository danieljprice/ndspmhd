!!-----------------------------------------------------------------------
!! Using the link list setup, compiles the neighbour list for the
!! current cell (this list is common to all particles in the cell)
!!
!! the list is returned in 'listneigh' (length nneigh)
!! the neigbouring cells to the current cell are returned in 'neighcell'
!!
!! This version for partial lists of particles (finds all neighbours)
!!-----------------------------------------------------------------------

SUBROUTINE get_neighbour_list_partial(icell,neighcell,listneigh,nneigh)     
 USE dimen_mhd 
 USE debug
 USE loguns
 USE linklist
 USE options
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: icell
 INTEGER, DIMENSION(3**ndim), INTENT(OUT) :: neighcell
 INTEGER, DIMENSION(nlistdim), INTENT(OUT) :: listneigh
 INTEGER, INTENT(OUT) :: nneigh
 INTEGER :: nneighcell,ipart
 INTEGER :: i,j,k

!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'

 neighcell(1) = icell		! first neighbouring cell (itself)
 neighcell(2) = icell + 1	! cell to the right
 neighcell(3) = icell - 1
 
 IF (ndim.GE.2) THEN
    neighcell(4) = icell + ncellsx(1)		! above
    neighcell(5) = icell + ncellsx(1) + 1	! above, one to right
    neighcell(6) = icell + ncellsx(1) - 1	! above, one to left
    neighcell(7) = icell - ncellsx(1)		! below
    neighcell(8) = icell - ncellsx(1) + 1	! below, one to right
    neighcell(9) = icell - ncellsx(1) - 1	! below, one to left
    IF (ndim.GE.3) THEN
        neighcell(10) = icell + ncellsx(1)*ncellsx(2)	 ! next block
	neighcell(11) = icell + ncellsx(1)*ncellsx(2) + 1 ! next block,right
	neighcell(12) = icell + ncellsx(1)*ncellsx(2) - 1 ! next block,left	
	neighcell(13) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) !  " ", above
	neighcell(14) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) + 1 ! " ", above,right
	neighcell(15) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) - 1 ! " ", above,left
    ENDIF
 ENDIF
 nneighcell = 3**ndim        ! number of cells to interact with

 IF (MOD(icell,ncellsx(1)).EQ.1) THEN	! cells next to xmin
    nneighcell = nneighcell - ndim
    IF (ndim.GE.2) THEN	! overwrite cells to left
       neighcell(3) = icell - ncellsx(1)
       neighcell(6) = icell - ncellsx(1) + 1
    ENDIF
 ELSEIF (MOD(icell,ncellsx(1)).EQ.0) THEN ! cells next to xmax
    nneighcell = nneighcell - ndim	! for rightmost cells, no cells to right
    neighcell(2) = icell - 1	!  one to left
    IF (ndim.GE.2) THEN
       neighcell(3) = icell - ncellsx(1)
       neighcell(5) = icell - ncellsx(1) - 1
    ENDIF
 ENDIF
! IF ((icell.EQ.1).AND.(ibound.GT.1)) THEN	! if first cell and 
!    nneighcell = 3				! ghosts set, they are in
!    neighcell(3) = icell - 1			! cell to left (cell 0)  		
! ELSEIF (icell.EQ.ncells) THEN
!    IF (ibound.GT.1) THEN	! if ghosts set, they are in ncells+1
!       nneighcell = 2		! last cell has itself and ghost cell
!    ELSE
!       nneighcell = 1		! for last cell, just itself
!    ENDIF
! ENDIF
	 
! PRINT*,'cell ',icell,' first particle = ',ifirstincell(icell)
!
!--construct list of neighbours for particles in the current cell
!	 
 j=0
 DO k = 1,nneighcell		! construct list of neighbours
    ipart = ifirstincell(neighcell(k))
!   PRINT*,'neighbouring cell = ',neighcell(k)
    IF (ipart.NE.-1) THEN
       DO WHILE (ll(ipart).NE.-1)	       
          j = j + 1
          IF (j.GT.nlistdim) THEN	! should reallocate array here
             WRITE(iprint,*) 'getneigh: # neighbours > array size'
	     CALL quit
          ENDIF
          listneigh(j) = ipart
          ipart = ll(ipart)
!         print*,'listneigh ',j,'= ',listneigh(j)
      ENDDO
      j = j + 1
      listneigh(j) = ipart
    ENDIF    
!   print*,'listneigh ',j,'= ',listneigh(j)	    
 ENDDO
 nneigh = j
!	 print*,'number of neighbours = ',nneigh
!	 read*	

 RETURN
END SUBROUTINE get_neighbour_list_partial
