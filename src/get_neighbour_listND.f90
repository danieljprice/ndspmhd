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
 INTEGER, DIMENSION(*), INTENT(OUT) :: listneigh
 INTEGER, INTENT(OUT) :: nneigh
 INTEGER, DIMENSION(3**ndim), INTENT(OUT) :: neighcell
 INTEGER :: nneighcell,ipart,ncellsxy
 INTEGER :: i,j,k
 LOGICAL :: debugging
 LOGICAL :: leftmost,rightmost,bottomrow,toprow,endblock
!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'
 debugging = .false.
 IF (idebug(1:4).EQ.'neig') debugging = .true.
!
!--work out whether cell is near a boundary
!
 leftmost = .false.
 rightmost = .false.
 toprow = .false.
 bottomrow = .false.
 endblock = .false.

 IF (MOD(icell,ncellsx(1)).EQ.1) leftmost = .true.
 IF (MOD(icell,ncellsx(1)).EQ.0) rightmost = .true.
 IF (ndim.GE.2) THEN
    ncellsxy = ncellsx(2)*ncellsx(1)
    IF (MOD(icell-1,ncellsxy).GT.(ncellsxy-ncellsx(1)-1)) toprow = .true.
    IF (MOD(icell-1,ncellsxy).LE.ncellsx(1)-1) bottomrow = .true.    
    IF (ndim.GE.3 .AND. (ncells-icell.LT.ncellsxy)) endblock = .true.    
 ENDIF
!
!--then work out which cells current cell should interact with
!
 nneighcell = 1
 neighcell(1) = icell ! always interacts with itself
 
 IF (.NOT.rightmost) THEN   ! cell to the right
    nneighcell = nneighcell + 1
    neighcell(nneighcell) = icell + 1
 ENDIF
 IF (ndim.GE.2 .AND..NOT.toprow) THEN ! above  
    IF (.NOT.leftmost) THEN   ! above, left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) - 1
    ENDIF
    nneighcell = nneighcell + 1   ! above
    neighcell(nneighcell) = icell + ncellsx(1)
    IF (.NOT.rightmost) THEN   ! above, right
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) + 1
    ENDIF
 ENDIF
 IF (ndim.EQ.3 .AND..NOT.endblock) THEN  ! next block
    IF (.NOT.leftmost) THEN    ! next block, left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy - 1    
    ENDIF   
    nneighcell = nneighcell + 1   ! next block
    neighcell(nneighcell) = icell + ncellsxy
    IF (.NOT.rightmost) THEN    ! next block, right  
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy + 1
    ENDIF
    IF (.NOT.toprow) THEN  ! next block, above
       IF (.NOT.leftmost) THEN   ! next block, above left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) - 1
       ENDIF
       nneighcell = nneighcell + 1  ! next block, above 
       neighcell(nneighcell) = icell + ncellsxy + ncellsx(1)
       IF (.NOT.rightmost) THEN   ! next block, above right
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) + 1
       ENDIF   
    ENDIF
    IF (.NOT.bottomrow) THEN  ! next block, below
       IF (.NOT.leftmost) THEN   ! next block, below left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) - 1
       ENDIF
       nneighcell = nneighcell + 1  ! next block, below 
       neighcell(nneighcell) = icell + ncellsxy - ncellsx(1)
       IF (.NOT.rightmost) THEN   ! next block, below right
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) + 1
       ENDIF   
    ENDIF
 ENDIF

 IF (debugging) THEN
    PRINT*,' current cell = ',icell,' first particle = ',ifirstincell(icell)
    IF (rightmost) PRINT*,'rightmost'
    IF (leftmost) PRINT*,'leftmost'
    IF (toprow) PRINT*,'toprow'
    IF (bottomrow) PRINT*,'bottomrow'
    IF (endblock) PRINT*,'endblock'
    PRINT*,' number of neighbouring cells = ',nneighcell
    PRINT*, (neighcell(i),i=1,nneighcell)
    READ*
 ENDIF  
!
!--construct list of neighbours for particles in the current cell
!  using the link lists from the neighbouring cells
!  
 j=0
 DO k = 1,nneighcell  ! construct list of neighbours
    IF (neighcell(k).GT.ncells) THEN
      PRINT*,' k, neighcell = ',k,neighcell(1:nneighcell),ncells
      STOP 'get_neighbour_list: error: cell > ncells'      
    ENDIF   
    ipart = ifirstincell(neighcell(k))
!   PRINT*,'neighbouring cell = ',neighcell(k)
    IF (ipart.NE.-1) THEN
       DO WHILE (ll(ipart).NE.-1)        
          j = j + 1
          IF (j.GT.nlistdim) THEN ! should reallocate array here
             WRITE(iprint,*) 'getneigh: # neighbours > array size:',nlistdim
             WRITE(iprint,*) 'ipart = ',ipart
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
!  read* 

 RETURN
END SUBROUTINE get_neighbour_list
