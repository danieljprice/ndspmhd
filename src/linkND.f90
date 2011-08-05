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

SUBROUTINE link
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE options
 USE part
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,icell
 INTEGER, DIMENSION(ndim) :: icellx
 REAL, DIMENSION(ndim) :: xminpart,xmaxpart
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine link'
!
!--use maximum value of h for cell size
! (this is already calculated in set_ghosts if using ghost particles)
 
 IF (ibound.LE.1) hhmax = MAXVAL(hh(1:npart))
 
 dxcell = radkern*hhmax		! size of link list cell (=compact support region)
 IF (dxcell.LE.0) THEN
    WRITE(iprint,*) 'Error: link: max h <=0 :',hhmax     
    CALL quit
 ENDIF
!
!--find max/min of particle distribution or use boundaries already set
!  
 xminpart(:) = xmin(:)
 xmaxpart(:) = xmax(:)
 
 IF (ibound.LE.1) THEN
!
!--find max/min of particle distribution, including ghost particles
!  (add/subtract 0.00001 to make sure all particles are within these limits)

    DO j=1,ndim
       xminpart(j) = MINVAL(x(j,1:ntotal)) - 0.00001
       xmaxpart(j) = MAXVAL(x(j,1:ntotal)) + 0.00001
    ENDDO
 ELSE
    xminpart(:) = xminpart(:) - dxcell - 0.00001    ! change boundaries to include ghosts
    xmaxpart(:) = xmaxpart(:) + dxcell + 0.00001
 ENDIF

 ncellsx(:) = int((xmaxpart(:)-xminpart(:))/dxcell) + 1
 
 IF (ANY(ncellsx .EQ. 0)) THEN
    WRITE(iprint,*) 'Error: link: number of cells=0:'
    WRITE(iprint,*) 'xmin, xmax, dxcell, hhmax = ',xmin,xmax,dxcell,hhmax
    WRITE(iprint,*) 'Max h particle ',MAXLOC(hh)
    CALL quit
 ENDIF 
! IF (ibound.GE.1) THEN	! if there are any boundaries
!    dxcell = (xmaxpart(1)-xminpart(1))/ncellsx(1)	! adjust so that xmax-xmin is exact division of 2h
! ELSE
!    ncellsx(:) = ncellsx(:) + 1
! ENDIF
 
 ncells = PRODUCT(ncellsx)

! WRITE(iprint,*) ' ncells x,y,z, hhmax = ',ncells,ncellsx,hhmax,dxcell

!
!--when doing calculation, in 1/2/3D 
!  don't do last cell/row of cells/face of cells as these are purely ghosts
!  really should make a list of cells to do in 3D as should skip the top row
!  of every block. At the moment this is just slightly inefficient.
!
 IF (ibound.LE.1) THEN	! BUG FIXED HERE (should be ibound.LE.1)
    ncellsloop = ncells
 ELSE
    IF (ndim.EQ.1) ncellsloop = ncells - 1
    IF (ndim.EQ.2) ncellsloop = ncells-ncellsx(1)
    IF (ndim.EQ.3) ncellsloop = ncells-ncellsx(1)*ncellsx(2) ! inefficient in 3D
 ENDIF
! print*,' ncellsloop = ',ncellsloop
! read*

!
! allocate array memory now we know the number of cells
! 
 DEALLOCATE( ifirstincell )
 ALLOCATE ( ifirstincell(ncells) )
 
! IF (ncells.GT.maxcells) PRINT*,'Error: link: ncells>array limits'
 DO i=1,ncells
    ifirstincell(i) = -1 		! set all head of chains to -1
 ENDDO

!
!--work out which cell the particle is in
!  
 DO i=1,ntotal		! including ghosts
    ll(i) =  -1
    icellx(:) = int((x(:,i)-xminpart(:))/dxcell) + 1 
    IF ( ANY(icellx < 0).OR.ANY(icellx > ncellsx) ) THEN
       WRITE(iprint,*) 'link: particle crossed boundary ',i,x(:,i),icellx,ncellsx
       CALL quit
    ENDIF
	 
    icell = icellx(1) 
    IF (ndim.GE.2) THEN
       j = 2
       icell = icell + (icellx(j)-1)*ncellsx(j-1)
       IF (ndim.GE.3) THEN
          j = 3
          icell = icell + (icellx(j)-1)*ncellsx(j-2)*ncellsx(j-1)
       ENDIF
    ENDIF

    ll(i) = ifirstincell(icell)		! link to previous start of chain
    ifirstincell(icell) = i		! set head of chain to current particle	 
    iamincell(i) = icell 		! save which cell particle is in 
!    PRINT*,' particle ',i,' x =',x(:,i),' in cell ',icellx,' = ',icell
!    read*
 ENDDO
!
!--debugging
!
 IF (idebug(1:4).EQ.'link') THEN
    print*,'xmin,xmax,dxcell,hhmax,ncells',xmin,xmax,dxcell,hhmax,ncells
    DO i=1,ncells
       j = ifirstincell(i)
       print*,'cell ',i,' ifirstincell = ',j
       DO WHILE (j.NE.-1 .AND. ll(j).NE.-1)
          j = ll(j)
          print*,' next =',j
       ENDDO
	 read*
    ENDDO     
 ENDIF    
END SUBROUTINE link
      
