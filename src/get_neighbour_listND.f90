!!-----------------------------------------------------------------------
!! Using the link list setup, compiles the neighbour list for the
!! current cell (this list is common to all particles in the cell)
!!
!! the list is returned in 'listneigh' (length nneigh)
!! the neigbouring cells to the current cell are returned in 'neighcell'
!!
!!-----------------------------------------------------------------------

SUBROUTINE get_neighbour_list(icell,neighcell,listneigh,nneigh)     
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
 INTEGER, DIMENSION(3*2**(ndim-1) - 1), INTENT(OUT) :: neighcell
 INTEGER, DIMENSION(*), INTENT(OUT) :: listneigh
 INTEGER, INTENT(OUT) :: nneigh
 INTEGER :: nneighcell,ipart
 INTEGER :: i,j,k

!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'

 neighcell(1) = icell		! first neighbouring cell (itself)
 neighcell(2) = icell + 1	! cell to the right
 
 IF (ndim.GE.2) THEN
    neighcell(3) = icell + ncellsx(1)		! above
    neighcell(4) = icell + ncellsx(1) + 1	! above, one to right
    neighcell(5) = icell + ncellsx(1) - 1	! above, one to left
    IF (ndim.GE.3) THEN
        neighcell(6) = icell + ncellsx(1)*ncellsx(2)	 ! next block
	neighcell(7) = icell + ncellsx(1)*ncellsx(2) + 1 ! next block,right
	neighcell(8) = icell + ncellsx(1)*ncellsx(2) - 1 ! next block,left	
	neighcell(9) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) !  " ", above
	neighcell(10) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) + 1 ! " ", above,right
	neighcell(11) = icell + ncellsx(1)*ncellsx(2) + ncellsx(1) - 1 ! " ", above,left
    ENDIF
 ENDIF
 nneighcell = 3*2**(ndim-1) - 1        ! number of cells to interact with

 IF (ndim.GE.2 .AND. MOD(icell,ncellsx(1)).EQ.1) THEN
    nneighcell = nneighcell - 1	! for leftmost cells, no cells above left
 ELSEIF (MOD(icell,ncellsx(1)).EQ.0) THEN
    nneighcell = nneighcell - ndim	! for rightmost cells, no cells to right
    neighcell(2) = icell + ncellsx(1) - 1	! above, one to left
 ENDIF
 
 IF (ndim.GE.2 .AND. icell.GE.ncells-ncellsx(1)) THEN
!!--for top cells and no boundaries, no cells above
    nneighcell = 2
    IF (icell.EQ.ncells) nneighcell = 1
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
    IF (neighcell(k).GT.ncells) THEN
      PRINT*,' k, neighcell = ',k,neighcell(1:nneighcell),ncells
      STOP 'error: cell > ncells'      
    ENDIF   
    ipart = ifirstincell(neighcell(k))
!   PRINT*,'neighbouring cell = ',neighcell(k)
    IF (ipart.NE.-1) THEN
       DO WHILE (ll(ipart).NE.-1)	       
          j = j + 1
          IF (j.GT.nlistdim) THEN	! should reallocate array here
             WRITE(iprint,*) 'getneigh: # neighbours > array size:',nlistdim
	     WRITE(iprint,*) 'ipart = ',ipart
	     WRITE(iprint,*) listneigh(1:nlistdim)
	     CALL quit
          ENDIF
          listneigh(j) = ipart
          ipart = ll(ipart)
!         print*,'listneigh ',j,'= ',listneigh(j)
      ENDDO
      j = j + 1
      listneigh(j) = ipart
    ENDIF    
!    IF (j.GT.0) THEN 
!       print*,'listneigh ',j,'= ',listneigh(j)	    
!    ELSE
!       print*,' no neighbours'
!    ENDIF   
 ENDDO
 nneigh = j
! print*,'number of neighbours = ',nneigh
!	 read*	

 RETURN
END SUBROUTINE get_neighbour_list
