!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab (h_a)
!!
!! Also computes the variable smoothing length terms sum_b m_b dWdh_b dhdrho_b
!!
!! This version computes the density on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

SUBROUTINE density
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
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,icellloop,ipart,iprev,ncell,nneigh
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh ! neighbour list
 INTEGER :: idone
 INTEGER, DIMENSION(3**ndim) :: neighcell
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2,rho1i
 REAL :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwab,hfacwabi,hfacwabj
 REAL, DIMENSION(ndim) :: dx
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj,weight
 REAL :: grkern,grkerni,grkernj
!
!  (grad h terms)
!
 REAL :: dwdhi,dwdhj,dhdrhoi,dhdrhoj
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine density'  
!
!--initialise quantities
!
 nlistdim = ntotal
 ALLOCATE( listneigh(nlistdim) )	! max size of neighbour list

 DO i=1,npart
    rho(i) = 0.
    gradh(i) = 0.
 ENDDO
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
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       dhdrhoi = -hi*dndim	! divide by  should use rho(i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
	  j = listneigh(n)
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
	  
	  rij2 = DOT_PRODUCT(dx,dx)
	  rij = SQRT(rij2)
	  q2 = rij2/h2
	  q2i = rij2/hi2
	  q2j = rij2/hj2	
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!	
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
	 IF (((q2.LT.radkern2).OR.(q2i.LT.radkern2).OR.(q2j.LT.radkern2))  &
	     .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!	
!--interpolate from kernel table		
!  (use either average h or average kernel gradient)
!
!	  PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
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
!
!--derivative w.r.t. h for grad h correction terms (and dhdrho)
!	       
	       dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
	       dwdhj = -rij*grkernj*hj1 - ndim*wabj*hj1
	       dhdrhoj = -hj*dndim	! see above		   				  
	    ENDIF
!
!--calculate density
!
	    weight = 1.0
	    IF (j.EQ.i) weight = 0.5
	    IF (ikernav.EQ.3) THEN
	       rho(i) = rho(i) + pmass(j)*wabi*weight
	       rho(j) = rho(j) + pmass(i)*wabj*weight
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

	       gradh(i) = gradh(i) + dhdrhoi*pmass(j)*weight*dwdhi
	       gradh(j) = gradh(j) + dhdrhoj*pmass(i)*weight*dwdhj
		  
	    ELSE
	       rho(i) = rho(i) + pmass(j)*wab*weight		
	       rho(j) = rho(j) + pmass(i)*wab*weight	
	    ENDIF
!	 ELSE
!	    PRINT*,' r/h > 2 '      
	 
	 ENDIF
	    
       ENDDO loop_over_neighbours

       iprev = i
       IF (iprev.NE.-1) i = ll(i)		! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells

 RETURN
END SUBROUTINE density
      
