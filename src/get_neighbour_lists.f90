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

  subroutine get_neighbour_list(icell,listneigh,nneigh)     
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
    integer, dimension(3**ndim) :: neighcell
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
!--if cell is empty, neighbouring cells are irrelevant, so return
!
    if (ifirstincell(icell).le.0) then
       if (debugging) write(iprint,*) 'icell:',icell,' cell empty: no neighbours'
       listneigh = 0
       nneigh = 0
       return
    endif
!
!--work out which cells current cell should interact with
!  same template is used for all cells as a block of empty cells surrounds those with particles
!  (ie. so there is always a cell to the left, above or whatever)
!
    nneighcell = 2
    neighcell(1) = icell ! always interacts with itself
    neighcell(2) = icell + 1
    if (ndim.ge.2) then
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) - 1   ! above, left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1)       ! above
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) + 1   ! above, right
    endif
    if (ndim.ge.3) then  ! next block
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy - 1    ! next block, left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy        ! next block
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy + 1    ! next block, right
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) - 1 ! next block, above left
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy + ncellsx(1)     ! next block, above
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) + 1 ! next block, above right
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) - 1 ! next block, below left
       nneighcell = nneighcell + 1  ! next block, below 
       neighcell(nneighcell) = icell + ncellsxy - ncellsx(1)     ! next block, below
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) + 1 ! next block, below right
    endif

    if (debugging) then
       print*,' current cell = ',icell,' first particle = ',ifirstincell(icell)
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

  subroutine get_neighbour_list_partial(icell,listneigh,nneigh)     
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
    integer, dimension(3**ndim) :: neighcell
    integer :: nneighcell,ipart,ncellsxy
    integer :: i,j,k, n2d, n3d,ineighcell
    logical :: debugging
    logical :: leftmost,rightmost,bottomrow,toprow,endblock,firstblock
!
!--allow for tracing flow (removed for speed)
!      
! IF (trace) WRITE(iprint,*) ' Entering subroutine get_neighbour_list'
    debugging = .false.
    if (idebug(1:4).EQ.'neig') debugging = .true.
!
!--if cell is empty, neighbouring cells are irrelevant, so return
!
    if (ifirstincell(icell).le.0) then
       if (debugging) write(iprint,*) 'icell:',icell,' cell empty: no neighbours'
       listneigh = 0
       nneigh = 0
       return
    endif
!
!--work out which cells current cell should interact with
!
    nneighcell = 3
    neighcell(1) = icell - 1  ! always interacts with itself
    neighcell(2) = icell
    neighcell(3) = icell + 1
    if (ndim.ge.2) then       ! add incrementally to stop compiler errors if ndim < 2
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) - 1
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1)
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell + ncellsx(1) + 1
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell - ncellsx(1) - 1
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell - ncellsx(1)
       nneighcell = nneighcell + 1
       neighcell(nneighcell) = icell - ncellsx(1) + 1
       if (ndim.ge.3) then
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1)
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + ncellsx(1) + 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy + 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1)
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell + ncellsxy - ncellsx(1) + 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy + ncellsx(1) - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy + ncellsx(1)
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy + ncellsx(1) + 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy + 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy - ncellsx(1) - 1
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy - ncellsx(1)
          nneighcell = nneighcell + 1
          neighcell(nneighcell) = icell - ncellsxy - ncellsx(1) + 1
       endif
    endif

    if (debugging) then
       print*,' current cell = ',icell,' first particle = ',ifirstincell(icell)
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
