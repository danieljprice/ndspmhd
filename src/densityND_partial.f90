!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab (h_a)
!!
!! This version computes the density only on a selected list of particles
!! given by the contents of the array ipartlist (enables iteration on 
!! unconverged particles only). It is therefore slightly slower
!! since some particle pairs may be done twice, once for each particle.
!!
!! Assumes rho_a is only a function of h_a
!!
!! This version must be used for individual particle timesteps
!!------------------------------------------------------------------------

SUBROUTINE density_partial(nlist,ipartlist)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE hterms
 USE kernel
 USE linklist
! USE options
 USE part
! USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,n
 INTEGER :: icell,icellloop,ipart,iprev,ncell,nneigh,index,index1
 INTEGER :: listneigh(nlistdim) ! up to 10% of particles in each cell
 INTEGER :: idone,icellprev
 INTEGER, DIMENSION(3**ndim) :: neighcell
 INTEGER, INTENT(IN) :: nlist
 INTEGER, INTENT(IN), DIMENSION(*) :: ipartlist
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2,rho1i
 REAL :: hi,hi1,hi2,hi3
 REAL :: hfacwabi,hfacgrkerni
 REAL, DIMENSION(ndim) :: dx
!
!  (kernel quantities)
!
 REAL :: q2i      
 REAL :: wabi,grkerni,weight      
!
!  (grad h terms)
!
 REAL :: dhdrhoi,dwdhi
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine density_partial'  
!
!--initialise quantities
!
 DO ipart=1,nlist
    i = ipartlist(ipart)
    rho(i) = 0.
    gradh(i) = 0.
 ENDDO
 icellprev = 0
!
!--Loop over all the particles in the density list
!
 loop_over_particles: DO ipart=1,nlist		! step through all cells

    i = ipartlist(ipart)
!
!--find cell of current particle
!
    icell = iamincell(i)
!    PRINT*,' particle ',i,' cell = ',icell
!
!--if different to previous cell used, get the list of neighbours for this cell
!  (common to all particles in the cell)
!
    IF (icell.NE.icellprev) THEN
       CALL get_neighbour_list_partial(icell,neighcell,listneigh,nneigh)
    ENDIF
    icellprev = icell

!   PRINT*,'Doing particle ',i,nneigh,' neighbours',pmass(i)
    hi = hh(i)
    hi2 = hi*hi
    hi3 = hi*hi2
    hi1 = 1./hi
       
    hfacwabi = hi1**ndim
    hfacgrkerni = hfacwabi*hi1
!
!--loop over current particle's neighbours
!
    loop_over_neighbours: DO n = 1,nneigh
       j = listneigh(n)
       dx(:) = x(:,i) - x(:,j)
!
!--calculate averages of smoothing length if using this averaging
!			 	  
       rij2 = DOT_PRODUCT(dx,dx)
       rij = SQRT(rij2)
       q2i = rij2/hi2
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2i),dx,hi
!	
!--do interaction if r/h < compact support size
!
       IF (q2i.LT.radkern2) THEN
!	
!--interpolate from kernel table (using hi)
!
          CALL interpolate_kernel(q2i,wabi,grkerni)
	  wabi = wabi*hfacwabi
	  grkerni = grkerni*hfacgrkerni
!
!--derivative w.r.t. h for grad h correction terms (and dhdrho)
!	       
          dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
	  dhdrhoi = -hi*dndim
!
!--calculate density
!
          rho(i) = rho(i) + pmass(j)*wabi
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

	  gradh(i) = gradh(i) + dhdrhoi*pmass(j)*dwdhi
	 
       ENDIF
	    
    ENDDO loop_over_neighbours

 ENDDO loop_over_particles

 RETURN
END SUBROUTINE density_partial
      
