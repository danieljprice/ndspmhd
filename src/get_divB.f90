!!------------------------------------------------------------------------
!! Computes an SPH estimate of div B
!! This version computes div B on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

SUBROUTINE get_divB(divBonrho,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE options
 USE part
 USE setup_params
 USE get_neighbour_lists
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ntot
 INTEGER :: i,j,n
 INTEGER :: icell,iprev,nneigh,nlistdim
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh ! neighbour list
 INTEGER :: idone
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2
 REAL :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwab,hfacwabi,hfacwabj
 REAL :: rho21i, rho21j, term, projdB
 REAL, DIMENSION(ndim) :: dx
 REAL, DIMENSION(ndimV) :: dr
 REAL, DIMENSION(ntot), INTENT(OUT) :: divBonrho
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj
 REAL :: grkern,grkerni,grkernj
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine get_divB, ntot=',ntot  
!
!--initialise quantities
!
 nlistdim = ntotal
 ALLOCATE( listneigh(nlistdim) )     ! max size of neighbour list

 DO i=1,ntot
    divBonrho(i) = 0.
 ENDDO
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop          ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)          ! loop over home cell particles

       !       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       rho21i = 1./rho(i)**2
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
          j = listneigh(n)
          IF ((j.NE.i).AND..NOT.(j.GT.npart .AND. i.GT.npart)) THEN
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             hj = hh(j)
             hj1 = 1./hj
             hj2 = hj*hj
!
!--calculate averages of smoothing length if using this averaging
!                
             hav = 0.5*(hi + hj)
             hav1 = 1./hav
             h2 = hav*hav
             hfacwab = hav1**ndim
             hfacwabj = hj1**ndim
             rho21j = 1./rho(j)**2
             
             rij2 = DOT_PRODUCT(dx,dx)
             rij = SQRT(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2     
             dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             if (ndimV.gt.ndim) dr(ndim+1:ndimV) = 0. 

             !          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             IF (((q2i.LT.radkern2).OR.(q2j.LT.radkern2))  &
                  .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
                IF (ikernav.EQ.1) THEN          
                   CALL interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hj1
                ELSE
                   !  (using hi)
                   CALL interpolate_kernel(q2i,wabi,grkerni)
                   wabi = wabi*hfacwabi
                   grkerni = grkerni*hfacwabi*hi1
                   !  (using hj)
                   CALL interpolate_kernel(q2j,wabj,grkernj)
                   wabj = wabj*hfacwabj
                   grkernj = grkernj*hfacwabj*hj1
                   !  (calculate average)            
                   wab = 0.5*(wabi + wabj)                  
                   grkern = 0.5*(grkerni + grkernj)            
                ENDIF
!
!--calculate div B
!
                term = DOT_PRODUCT((Bfield(:,i)*rho21i + Bfield(:,j)*rho21j),dr)
                projdB = DOT_PRODUCT(Bfield(:,i)-Bfield(:,j),dr)
                IF (ikernav.EQ.3) THEN
                   !               divBonrho(i) = divBonrho(i) - pmass(j)*projdB*grkerni
                   !               divBonrho(j) = divBonrho(j) - pmass(i)*projdB*grkernj            
                   
                   divBonrho(i) = divBonrho(i) + pmass(j)*term*grkerni
                   divBonrho(j) = divBonrho(j) - pmass(i)*term*grkernj            
                ELSE
                   !               divBonrho(i) = divBonrho(i) - pmass(j)*projdB*grkern
                   !               divBonrho(j) = divBonrho(j) - pmass(i)*projdB*grkern            
                   
                   divBonrho(i) = divBonrho(i) + pmass(j)*term*grkern     
                   divBonrho(j) = divBonrho(j) - pmass(i)*term*grkern     
                ENDIF
                !      ELSE
                !         PRINT*,' r/h > 2 '      
                
             ENDIF
          ENDIF! j .ne. i   
       ENDDO loop_over_neighbours
       
       iprev = i
       IF (iprev.NE.-1) i = ll(i)          ! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
    
 ENDDO loop_over_cells

! do i=1,ntotal
!    divBonrho(i) = divBonrho(i)/rho(i)**2
! enddo

 RETURN
END SUBROUTINE get_divB
      
