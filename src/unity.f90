!!------------------------------------------------------------------------
!! Computes the unity function by direct summation over the particles
!! ie. unity_a = sum_b m_b/rho_b W_ab (h_a)
!! and also its derivative.
!! 
!! Density must already be known (either from the continuity equation or
!! by direct summation)
!!
!! This version computes the density on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

SUBROUTINE calc_unity
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE hterms
 USE kernel
 USE linklist
 USE options
 USE part
 USE setup_params
 USE unityfunc
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,iprev,nneigh
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh ! neighbour list
 INTEGER :: idone
 INTEGER, DIMENSION(3**ndim) :: neighcell
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2
 REAL :: hi,hi1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwabi,hfacwabj,rho1i,rho1j
 REAL, DIMENSION(ndim) :: dx, dr
!
!  (kernel quantities)
!
 REAL :: q2i,q2j      
 REAL :: wabi,wabj,weight
 REAL :: grkerni,grkernj
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine unity'  
!
!--initialise quantities
!
 nlistdim = ntotal
 ALLOCATE( listneigh(nlistdim) )       ! max size of neighbour list

 unity = 0.
 gradunity = 0.
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop              ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1       ! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)              ! loop over home cell particles

!       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       rho1i = 1./rho(i)
       hfacwabi = hi1**ndim
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
         j = listneigh(n)
         dx(:) = x(:,i) - x(:,j)
         hj = hh(j)
         hj1 = 1./hj
         hj2 = hj*hj
	 rho1j = 1./rho(j)
         hfacwabj = hj1**ndim
         
         rij2 = DOT_PRODUCT(dx,dx)
         rij = SQRT(rij2)
	 if (rij.gt.1.e-8) then
	    dr(:) = dx(:)/rij
	 else
	    dr = 0.
	 endif
         q2i = rij2/hi2
         q2j = rij2/hj2       
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2i),dx,hi,hj
!       
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
        IF ((q2i.LT.radkern2).OR.(q2j.LT.radkern2)  &
            .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!       
!--interpolate from kernel table              
!  (use either average h or average kernel gradient)
!
!  (using hi)
           CALL interpolate_kernel(q2i,wabi,grkerni)
           wabi = wabi*hfacwabi
           grkerni = grkerni*hfacwabi*hi1
!  (using hj)
           CALL interpolate_kernel(q2j,wabj,grkernj)
           wabj = wabj*hfacwabj
           grkernj = grkernj*hfacwabj*hj1
!
!--calculate unity function density
!
           weight = 1.0
           IF (j.EQ.i) weight = 0.5
           unity(i) = unity(i) + pmass(j)*rho1j*wabi*weight
           unity(j) = unity(j) + pmass(i)*rho1i*wabj*weight
	   gradunity(:,i) = gradunity(:,i) + pmass(j)*rho1j*dr(:)*grkerni
	   gradunity(:,j) = gradunity(:,j) - pmass(i)*rho1i*dr(:)*grkernj
	   
!        ELSE
!           PRINT*,' r/h > 2 '      
        
         ENDIF
           
       ENDDO loop_over_neighbours

       iprev = i
       IF (iprev.NE.-1) i = ll(i)              ! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells

! PRINT*,'normalising density'
 IF (icty.eq.0) THEN
    DO i=1,npart
!!    print*,i,rho(i),unity(i),gradunity(:,i),rho(i)/unity(i)
       rho(i) = rho(i)/unity(i)
    ENDDO
 ENDIF
 
 DO i=npart+1,ntotal
    j = ireal(i)
    if (icty.eq.0) rho(i) = rho(j)
    unity(i) = unity(j)
    gradunity(:,i) = gradunity(:,j)
 ENDDO
!! read*

 RETURN
END SUBROUTINE calc_unity
      
