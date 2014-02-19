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

!
! utility to check whether the neighbour finding is working or not
! checks the neighbours 
!
subroutine check_neighbourlist
 use dimen_mhd, only:ndim
 use kernels, only:radkern2
 use part, only:x, hh, npart,ntotal
 use linklist
 use get_neighbour_lists
 implicit none
 integer :: i,j,n,ineigh
 integer :: icell,iprev,nneigh,icellprev
 integer, allocatable, dimension(:) :: listneigh
 integer, allocatable, dimension(:,:) :: neighbour_list
 integer :: idone, nlistdim
 real :: rij2,hi,q2i,q2j
 real, dimension(ndim) :: dx
 logical :: iok, partial

 partial = .true.

 print*,'*** checking neighbour search ***'
 if (partial) print*,'*** with get_neighbours_partial ***'
!
!--setup link list
!
 call set_linklist
!
!--perform neighbour interactions as usual
!
!
!--initialise quantities
!
 nlistdim = ntotal
 allocate( listneigh(nlistdim),numneigh(ntotal) )        ! max size of neighbour list
 allocate( neighbour_list(ntotal,nlistdim) )
 numneigh(:) = 0
 neighbour_list = 0
 
 if (.not.partial) then
!
!--Loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop                ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1        ! note density summation includes current particle
    if (i.NE.-1) iprev = i

    loop_over_cell_particles: do while (i.NE.-1)                ! loop over home cell particles

!       print*,'Doing particle ',i,nneigh,' neighbours',pmass(i)
       idone = idone + 1
!
!--loop over neighbours as usual
! 
       hi = hh(i)
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          dx(:) = x(:,i) - x(:,j)
          rij2 = DOT_PRODUCT(dx,dx) 
          q2i = rij2/hi**2
          q2j = rij2/hh(j)**2        
          if ((q2i.LT.radkern2).OR.(q2j.LT.radkern2)) then  
             numneigh(i) = numneigh(i) + 1
             IF (i.NE.j) numneigh(j) = numneigh(j) + 1
             neighbour_list(i,numneigh(i)) = j
             neighbour_list(j,numneigh(j)) = i
          endif
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.NE.-1) i = ll(i)                ! possibly should be only IF (iprev.NE.-1)

    enddo loop_over_cell_particles
 enddo loop_over_cells
 
 else
!
!--do loop using get_neighbour_list_partial
!
    icellprev = 0
!
!--Loop over all the particles in the density list
!
    over_particles: do i=1,npart            ! step through all cells
!
!--find cell of current particle
!
       icell = iamincell(i)
!    PRINT*,' particle ',i,' cell = ',icell
!
!--if different to previous cell used, get the list of neighbours for this cell
!  (common to all particles in the cell)
!
       if (icell.NE.icellprev) then
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       endif
       icellprev = icell

       hi = hh(i)
       over_neighbours: do n = 1,nneigh
          j = listneigh(n)
          dx(:) = x(:,i) - x(:,j)
          rij2 = dot_product(dx,dx)
          q2i = rij2/hi**2
          if (q2i.LT.radkern2) then
             numneigh(i) = numneigh(i) + 1
             neighbour_list(i,numneigh(i)) = j
          endif
          
       enddo over_neighbours

    enddo over_particles
 endif
!
!--loop over all particles again
!
 do i=1,npart
    print*,'checking particle ',i,' number of neighbours = ',numneigh(i), 'cell    = ',iamincell(i)
    iok = .true.
!
!--now for this particle loop over every other particle and check 
!  each one individually to see if it should be a neighbour or not
! 
    hi = hh(i)
    ineigh = 0
    do j=1,ntotal
       dx(:) = x(:,i)-x(:,j)
       rij2 = DOT_PRODUCT(dx,dx)
       q2i = rij2/hi**2
       q2j = rij2/hh(j)**2
       if (q2i.LT.radkern2 .OR. q2j.LT.radkern2) then
          ineigh = ineigh + 1
          if (any(neighbour_list(i,:).EQ.j)) then
!             print*,j,' ok',iamincell(j)
          else
             print*,' neighbour ',j,' not found, x=',x(:,j),iamincell(j)
             iok = .false.
          endif
       endif
    enddo
    if (ineigh.NE.numneigh(i)) then
       print*,'error: should have ',ineigh,', neighbours, found ',numneigh(i)
!       print*,neighbour_list(i,1:numneigh(i))
    endif   
!    read*
    if (.not.iok) read*
 enddo
    
 return
end subroutine check_neighbourlist
