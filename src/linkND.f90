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

!!--------------------------------------------------------------------------
!! Creates link list to find particle neighbours (ND)
!!
!! Link list cells are numbered from xmin to xmax, then ymin to ymax
!! then zmin to zmax. All the particles are binned into link list 
!! cells, including ghosts.
!!
!! During the computation, each cell finds 2**ndim neighbouring cells
!! e.g. 1D, cell finds the cell to the right
!!          interactions between particles in both cells are then computed
!!          then move to next cell, calculate cell to the right.
!!      2D, cell finds 5 neighbouring cells
!!      3D, cell finds 11 neighbouring cells
!!
!! We only need to loop to ncells-1 in each direction, since each cell
!! finds the cell to the right as its neighbour and all the interactions
!! are computed.
!!
!! Changes log:
!! 17/10/03 - Bug fix for 1D fixed particle boundaries
!!--------------------------------------------------------------------------

subroutine set_linklist
! use dimen_mhd
 use debug, only:trace,idebug
 use loguns, only:iprint
 
 use bound, only:xmin,xmax,hhmax
 use kernels, only:radkern
 use linklist
 use options, only:ibound
 use part, only:x,hh,npart,ntotal
!
!--define local variables
!
 implicit none
 integer :: i,j,icell
 integer, dimension(ndim) :: icellx
 real, dimension(ndim) :: xminpart,xmaxpart
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine link'
!
!--use maximum value of h for cell size
! (this is already calculated in set_ghosts if using ghost particles)
 
 if (all(ibound.le.1)) hhmax = maxval(hh(1:npart))
 if (hhmax.gt.1.e10) write(iprint,*) 'ERROR: hmax > 10**10 on ',maxloc(hh(1:npart))
 dxcell = radkern*hhmax            ! size of link list cell (=compact support region)
 if (dxcell.le.0) then
    write(iprint,*) 'error: link: max h <=0 :',hhmax,maxloc(hh(1:npart))    
    call quit
 endif
!
!--find max/min of particle distribution, including ghost particles
!  (add/subtract 0.00001 to make sure all particles are within these limits)
!
 do j=1,ndim
    xminpart(j) = minval(x(j,1:ntotal)) - 0.00001
    xmaxpart(j) = maxval(x(j,1:ntotal)) + 0.00001
 enddo
!
!--add one (empty) cell to each end in each dimension
!
 xminpart(:) = xminpart(:) - dxcell - 0.00001
 xmaxpart(:) = xmaxpart(:) + dxcell + 0.00001
!
!--now work out total number of cells
!
 ncellsx(:) = int((xmaxpart(:)-xminpart(:))/dxcell) + 1
 if (any(ncellsx .eq. 0)) then
    write(iprint,*) 'error: link: number of cells=0:'
    write(iprint,*) 'xmin, xmax, dxcell, hhmax = ',xmin,xmax,dxcell,hhmax
    write(iprint,*) 'max h particle ',maxloc(hh)
    call quit
 endif
 
 ncells = product(ncellsx)

! write(iprint,*) ' ncells x,y,z, hhmax = ',ncells,ncellsx,hhmax,dxcell

 ncellsloop = ncells
!
! allocate array memory now we know the number of cells
! 
 if (allocated(ifirstincell)) deallocate( ifirstincell )
 allocate ( ifirstincell(ncells) )
 
 do i=1,ncells
    ifirstincell(i) = -1             ! set all head of chains to -1
 enddo

!
!--work out which cell the particle is in
!  
 do i=1,ntotal            ! including ghosts
    ll(i) =  -1
    icellx(:) = int((x(:,i)-xminpart(:))/dxcell) + 1 
    if ( any(icellx < 0).or.any(icellx > ncellsx) ) then
       write(iprint,*) 'link: particle crossed boundary ',i,x(:,i),icellx,ncellsx
       call quit
    endif
       
    icell = icellx(1) 
    if (ndim.ge.2) then
       j = 2
       icell = icell + (icellx(j)-1)*ncellsx(j-1)
       if (ndim.ge.3) then
          j = 3
          icell = icell + (icellx(j)-1)*ncellsx(j-2)*ncellsx(j-1)
       endif
    endif

    ll(i) = ifirstincell(icell)            ! link to previous start of chain
    ifirstincell(icell) = i            ! set head of chain to current particle       
    iamincell(i) = icell             ! save which cell particle is in 
    !if (i.eq.1) print*,' particle ',i,' x =',x(:,i),' in cell ',icellx,' = ',icell
 enddo
!
!--debugging
!
 if (idebug(1:4).eq.'link') then
    print*,'xmin,xmax,dxcell,hhmax,ncells',xmin,xmax,dxcell,hhmax,ncells
    do i=1,ncells
       j = ifirstincell(i)
       print*,'cell ',i,' ifirstincell = ',j
       do while (j.ne.-1 .and. ll(j).ne.-1)
          j = ll(j)
          print*,' next =',j
       enddo
       read*
    enddo     
 endif    
end subroutine set_linklist
      
