!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Smoothes out initial conditions using the SPH summation, that is the  !!
!!   quantity x_a on particle a is recalculated by summing over the       !!
!!   particle neighbours according to:                                    !!
!!                                                                        !!
!!   x_smooth(r_a) = sum_b m_b (x_b / rho_b) W(r-r_b,h_ab)                 !!
!!                                                                        !!
!!   where W is the SPH smoothing kernel 				  !!
!!									  !!
!!  Needs to be preceded by a call to linklist to setup neighbour lists   !!
!!  The density should also be previously known				  !!
!!------------------------------------------------------------------------!!

SUBROUTINE smooth_initial_conditions(xinput,xsmooth,nsize)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE kernel
 USE linklist
 USE part
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh ! neighbour list
 INTEGER :: idone,iprev,nsize,nneigh,icell
 INTEGER, DIMENSION(3*2**(ndim-1) - 1) :: neighcell
 REAL, DIMENSION(nsize), INTENT(IN) :: xinput
 REAL, DIMENSION(nsize), INTENT(OUT) :: xsmooth
 REAL, DIMENSION(ndim) :: dx
 REAL :: termi,termj,hav,hi,hav1,h21,hfacwab
 REAL :: q2,rij2,weight,wab,grkern
!
!--allow for tracing flow
!
 WRITE(iprint,*) 'Smoothing initial conditions...'
!
!--initialise quantities
!
 nlistdim = ntotal
 ALLOCATE( listneigh(nlistdim) )	! max size of neighbour list

 xsmooth(:) = 0.
!
!--check size of input array is consistent with rho
!
 IF (SIZE(xinput).NE.SIZE(rho)) THEN 
    WRITE(iprint,*)  &
      'Error: smooth: input array and rho have different dimensions'
    CALL quit   
 ENDIF   
!
!--Loop over all the link-list cells
!
 loop_over_cells: DO icell=1,ncellsloop		! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    CALL get_neighbour_list(icell,neighcell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1	! note density summation includes current particle
    IF (i.NE.-1) iprev = i

    loop_over_cell_particles: DO WHILE (i.NE.-1)		! loop over home cell particles

!       PRINT*,'Doing particle ',i,nneigh,' neighbours',pmass(i)
       idone = idone + 1
       hi = hh(i)
       termi = pmass(i)*xinput(i)/rho(i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
	  j = listneigh(n)
	  dx(:) = x(:,i) - x(:,j)
	  termj = pmass(j)*xinput(j)/rho(j)
!
!--calculate averages of smoothing lengths
!			 
	  hav = 0.5*(hi + hh(j))
	  hav1 = 1./hav
	  h21 = hav1*hav1
	  hfacwab = hav1**ndim
	  
	  rij2 = DOT_PRODUCT(dx,dx)
	  q2 = rij2*h21
!	
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
	 IF ((q2.LT.radkern2)  &
	     .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!	
!--interpolate from kernel table (use average h)
!
             CALL interpolate_kernel(q2,wab,grkern)
	     wab = wab*hfacwab
!
!--calculate smoothed quantity
!
	    weight = 1.0
	    IF (j.EQ.i) weight = 0.5

            xsmooth(i) = xsmooth(i) + termj*wab*weight
	    xsmooth(j) = xsmooth(j) + termi*wab*weight
	 
	 ENDIF
	    
       ENDDO loop_over_neighbours
	    
       iprev = i
       IF (iprev.NE.-1) i = ll(i)		! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells
 
! WHERE (rho > 0.) 
!    xsmooth(:) = xsmooth(:)/rho(:)
! END WHERE
!
!--update ghosts
! 
 IF (ntotal.GT.npart) THEN
    DO i=npart+1,ntotal
       j = ireal(i)
       xsmooth(i) = xsmooth(j)
    ENDDO
 ENDIF   

 RETURN
END SUBROUTINE smooth_initial_conditions
