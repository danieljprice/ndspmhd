module get_neighbour_lists
  implicit none

contains

!!-----------------------------------------------------------------------
!! Using the link list setup, compiles the neighbour list for the
!! current cell (this list is common to all particles in the cell)
!!
!! the list is returned in 'listneigh' (length nneigh)
!! the neigbouring cells to the current cell are returned in 'neighcell'
!!
!!-----------------------------------------------------------------------

  subroutine get_neighbour_list(icell,neighcell,listneigh,nneigh)     
    use dimen_mhd 
    use debug
    use loguns
    use linklist
    use options
!
!--define local variables
!
    implicit none
    integer, intent(IN) :: icell
    integer, dimension(:), intent(OUT) :: listneigh
    integer, intent(OUT) :: nneigh
    integer, dimension(3**ndim), intent(OUT) :: neighcell
    integer :: nneighcell,ipart,ncellsxy
    integer :: i,j,k
    logical :: debugging
    logical :: leftmost,rightmost,bottomrow,toprow,endblock
!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'
    debugging = .false.
    if (idebug(1:4).EQ.'neig') debugging = .true.
!
!--work out whether cell is near a boundary
!
    leftmost = .false.
    rightmost = .false.
    toprow = .false.
    bottomrow = .false.
    endblock = .false.

    if (mod(icell,ncellsx(1)).EQ.1) leftmost = .true.
    if (mod(icell,ncellsx(1)).EQ.0) rightmost = .true.
    if (ndim.GE.2) then
       ncellsxy = ncellsx(2)*ncellsx(1)
       if (mod(icell-1,ncellsxy).GT.(ncellsxy-ncellsx(1)-1)) toprow = .true.
       if (mod(icell-1,ncellsxy).LE.ncellsx(1)-1) bottomrow = .true.    
       if (ndim.GE.3 .AND. (ncells-icell.LT.ncellsxy)) endblock = .true.    
    endif
!
!--then work out which cells current cell should interact with
!
    nneighcell = 1
    neighcell(1) = icell ! always interacts with itself
 
    if (.not.rightmost) then   ! cell to the right
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + 1
    endif
    if (ndim.GE.2 .AND..not.toprow) then ! above  
       if (.not.leftmost) then   ! above, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsx(1) - 1
       endif
       nneighcell = nneighcell + 1   ! above
       neighcell(nneighcell) = icell + ncellsx(1)
       if (.not.rightmost) then   ! above, right
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsx(1) + 1
       endif
    endif
    if (ndim.EQ.3 .AND..not.endblock) then  ! next block
       if (.not.leftmost) then    ! next block, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - 1    
       endif
       nneighcell = nneighcell + 1   ! next block
       neighcell(nneighcell) = icell + ncellsxy
       if (.not.rightmost) then    ! next block, right  
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + 1
       endif
       if (.not.toprow) then  ! next block, above
          if (.not.leftmost) then   ! next block, above left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1  ! next block, above 
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1)
          if (.not.rightmost) then   ! next block, above right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) + 1
          endif
       endif
       if (.not.bottomrow) then  ! next block, below
          if (.not.leftmost) then   ! next block, below left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1  ! next block, below 
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1)
          if (.not.rightmost) then   ! next block, below right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) + 1
          endif
       endif
    endif

    if (debugging) then
       print*,' current cell = ',icell,' first particle = ',ifirstincell(icell)
       if (rightmost) print*,'rightmost'
       if (leftmost) print*,'leftmost'
       if (toprow) print*,'toprow'
       if (bottomrow) print*,'bottomrow'
       if (endblock) print*,'endblock'
       print*,' number of neighbouring cells = ',nneighcell
       print*, (neighcell(i),i=1,nneighcell)
       read*
    endif
!
!--construct list of neighbours for particles in the current cell
!  using the link lists from the neighbouring cells
!  
    j=0
    do k = 1,nneighcell  ! construct list of neighbours
       if (neighcell(k).GT.ncells) then
          print*,' k, neighcell = ',k,neighcell(1:nneighcell),ncells
          stop 'get_neighbour_list: error: cell > ncells'      
       endif
       ipart = ifirstincell(neighcell(k))
!   PRINT*,'neighbouring cell = ',neighcell(k)
       if (ipart.NE.-1) then
          do while (ll(ipart).NE.-1)        
             j = j + 1
             if (j.GT.size(listneigh)) then ! should reallocate array here
                write(iprint,*) 'getneigh: # neighbours > array size:',size(listneigh)
                write(iprint,*) 'ipart = ',ipart
                call quit
             endif
             listneigh(j) = ipart
             ipart = ll(ipart)
!         print*,'listneigh ',j,'= ',listneigh(j)
          enddo
          j = j + 1
          listneigh(j) = ipart
       endif
!    IF (j.GT.0) THEN 
!       print*,'listneigh ',j,'= ',listneigh(j)     
!    ELSE
!       print*,' no neighbours'
!    ENDIF   
    enddo
    nneigh = j
! print*,'number of neighbours = ',nneigh
!  read* 

    return
  end subroutine get_neighbour_list

!!-----------------------------------------------------------------------
!! Using the link list setup, compiles the neighbour list for the
!! current cell (this list is common to all particles in the cell)
!!
!! the list is returned in 'listneigh' (length nneigh)
!! the neigbouring cells to the current cell are returned in 'neighcell'
!!
!!-----------------------------------------------------------------------

  subroutine get_neighbour_list_partial(icell,neighcell,listneigh,nneigh)     
    use dimen_mhd 
    use debug
    use loguns
    use linklist
    use options
!
!--define local variables
!
    implicit none
    integer, intent(IN) :: icell
    integer, dimension(:), intent(OUT) :: listneigh
    integer, intent(OUT) :: nneigh
    integer, dimension(3**ndim), intent(OUT) :: neighcell
    integer :: nneighcell,ipart,ncellsxy
    integer :: i,j,k
    logical :: debugging
    logical :: leftmost,rightmost,bottomrow,toprow,endblock,firstblock
!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'
    debugging = .false.
    if (idebug(1:4).EQ.'neig') debugging = .true.
!
!--work out whether cell is near a boundary
!
    leftmost = .false.
    rightmost = .false.
    toprow = .false.
    bottomrow = .false.
    endblock = .false.
    firstblock = .false.

    if (mod(icell,ncellsx(1)).EQ.1) leftmost = .true.
    if (mod(icell,ncellsx(1)).EQ.0) rightmost = .true.
    if (ndim.GE.2) then
       ncellsxy = ncellsx(2)*ncellsx(1)
       if (mod(icell-1,ncellsxy).GT.(ncellsxy-ncellsx(1)-1)) toprow = .true.
       if (mod(icell-1,ncellsxy).LE.ncellsx(1)-1) bottomrow = .true.    
       if (ndim.GE.3 .AND. (ncells-icell.LT.ncellsxy)) endblock = .true.  
       if (ndim.GE.3 .AND. (icell.LE.ncellsxy)) firstblock = .true.  
    endif
!
!--then work out which cells current cell should interact with
!  (loop from previous block 
!
    n3d = 1
    n2d = 1
    if (ndim.ge.2) then
       n2d = 3
       ineighcell = ineighcell - ncellsx(1)
    endif
    if (ndim.eq.3) then
       n3d = 3
       ineighcell = ineighcell - ncellsxy
    endif
    ineighcell = icell - ncellsxy - ncellsx(1) - 2
    
    do k=1,n3d
       ineighcell = ineighcell + ncellsxy
       do j=1,n2d
          ineighcell = ineighcell + ncellsx(1)
          do i=1,3
             ineighcell = ineighcell + 1
             
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = ineighcell
          enddo
       enddo
    enddo

    nneighcell = 1
    neighcell(1) = icell - 1  ! always interacts with itself
    neighcell(2) = icell
    neighcell(3) = icell + 1
    if (ndim.ge.2) then
       neighcell(4) = icell + ncellsx(1) - 1
       neighcell(5) = icell + ncellsx(1)
       neighcell(6) = icell + ncellsx(1) + 1
       neighcell(7) = icell - ncellsx(1) - 1
       neighcell(8) = icell - ncellsx(1)
       neighcell(9) = icell - ncellsx(1) + 1
       if (ndim.ge.3) then
          neighcell(10) = icell + ncellsxy + ncellsx(1) - 1
          neighcell(11) = icell + ncellsxy + ncellsx(1)
          neighcell(12) = icell + ncellsxy + ncellsx(1) + 1
          neighcell(13) = icell + ncellsxy - 1
          neighcell(14) = icell + ncellsxy
          neighcell(15) = icell + ncellsxy + 1
          neighcell(16) = icell + ncellsxy - ncellsx(1) - 1
          neighcell(17) = icell + ncellsxy - ncellsx(1)
          neighcell(18) = icell + ncellsxy - ncellsx(1) + 1
          neighcell(19) = icell - ncellsxy + ncellsx(1) - 1
          neighcell(20) = icell - ncellsxy + ncellsx(1)
          neighcell(21) = icell - ncellsxy + ncellsx(1) + 1
          neighcell(22) = icell - ncellsxy - 1
          neighcell(23) = icell - ncellsxy
          neighcell(24) = icell - ncellsxy + 1
          neighcell(25) = icell - ncellsxy - ncellsx(1) - 1
          neighcell(26) = icell - ncellsxy - ncellsx(1)
          neighcell(27) = icell - ncellsxy - ncellsx(1) + 1
       endif
    endif
    
 
    if (.not.rightmost) then         ! cell to the right
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + 1
    endif
    if (.not.leftmost) then         ! cell to the left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell - 1
    endif
    if (ndim.GE.2 .AND..not.toprow) then   ! above      
       if (.not.leftmost) then         ! above, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsx(1) - 1
       endif
       nneighcell = nneighcell + 1         ! above
       neighcell(nneighcell) = icell + ncellsx(1)
       if (.not.rightmost) then         ! above, right
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsx(1) + 1
       endif
    elseif (ndim.GE.2 .AND..not.bottomrow) then
       if (.not.leftmost) then         ! below, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsx(1) - 1
       endif
       nneighcell = nneighcell + 1         ! below
       neighcell(nneighcell) = icell - ncellsx(1)
       if (.not.rightmost) then         ! below, right
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsx(1) + 1
       endif
    endif
    if (ndim.EQ.3 .AND..not.endblock) then    ! next block
       if (.not.leftmost) then          ! next block, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - 1    
       endif
       nneighcell = nneighcell + 1         ! next block
       neighcell(nneighcell) = icell + ncellsxy
       if (.not.rightmost) then          ! next block, right  
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + 1
       endif
       if (.not.toprow) then      ! next block, above
          if (.not.leftmost) then         ! next block, above left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1      ! next block, above   
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1)
          if (.not.rightmost) then         ! next block, above right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) + 1
          endif
       endif
       if (.not.bottomrow) then      ! next block, below
          if (.not.leftmost) then         ! next block, below left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1      ! next block, below   
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1)
          if (.not.rightmost) then         ! next block, below right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) + 1
          endif
       endif
       elseif (ndim.EQ.3 .AND..not.firstblock) then
       if (.not.leftmost) then          ! previous block, left
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy - 1    
       endif
       nneighcell = nneighcell + 1         ! previous block
       neighcell(nneighcell) = icell - ncellsxy
       if (.not.rightmost) then          ! previous block, right  
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy + 1
       endif
       if (.not.toprow) then      ! previous block, above
          if (.not.leftmost) then         ! previous block, above left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell - ncellsxy + ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1      ! previous block, above   
          neighcell(nneighcell) = icell - ncellsxy + ncellsx(1)
          if (.not.rightmost) then         ! previous block, above right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell - ncellsxy + ncellsx(1) + 1
          endif
       endif
       if (.not.bottomrow) then      ! previous block, below
          if (.not.leftmost) then         ! previous block, below left
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell - ncellsxy - ncellsx(1) - 1
          endif
          nneighcell = nneighcell + 1      ! previous block, below   
          neighcell(nneighcell) = icell - ncellsxy - ncellsx(1)
          if (.not.rightmost) then         ! previous block, below right
             nneighcell = nneighcell + 1
             neighcell(nneighcell) = icell - ncellsxy - ncellsx(1) + 1
          endif
       endif
    endif

    if (debugging) then
       print*,' current cell = ',icell,' first particle = ',ifirstincell(icell)
       if (rightmost) print*,'rightmost'
       if (leftmost) print*,'leftmost'
       if (toprow) print*,'toprow'
       if (bottomrow) print*,'bottomrow'
       if (endblock) print*,'endblock'
       print*,' number of neighbouring cells = ',nneighcell
       print*, (neighcell(i),i=1,nneighcell)
       read*
    endif
!
!--construct list of neighbours for particles in the current cell
!  using the link lists from the neighbouring cells
!    
    j=0
    do k = 1,nneighcell      ! construct list of neighbours
       if (neighcell(k).GT.ncells) then
          print*,' k, neighcell = ',k,neighcell(1:nneighcell),ncells
          stop 'get_neighbour_list: error: cell > ncells'      
       endif
       ipart = ifirstincell(neighcell(k))
!   PRINT*,'neighbouring cell = ',neighcell(k)
       if (ipart.NE.-1) then
          do while (ll(ipart).NE.-1)          
             j = j + 1
             if (j.GT.size(listneigh)) then   ! should reallocate array here
                write(iprint,*) 'getneigh: # neighbours > array size:',size(listneigh)
                write(iprint,*) 'ipart = ',ipart
                call quit
             endif
             listneigh(j) = ipart
             ipart = ll(ipart)
!         print*,'listneigh ',j,'= ',listneigh(j)
          enddo
          j = j + 1
          listneigh(j) = ipart
       endif
!    IF (j.GT.0) THEN 
!       print*,'listneigh ',j,'= ',listneigh(j)       
!    ELSE
!       print*,' no neighbours'
!    ENDIF   
    enddo
    nneigh = j
! print*,'number of neighbours = ',nneigh
!    read*   

    return
  end subroutine get_neighbour_list_partial

end module get_neighbour_lists
