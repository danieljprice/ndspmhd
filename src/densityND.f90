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
 INTEGER :: icell,iprev,nneigh
 INTEGER, ALLOCATABLE, DIMENSION(:) :: listneigh ! neighbour list
 INTEGER :: idone
 INTEGER, DIMENSION(3**ndim) :: neighcell
!
!  (particle properties - local copies)
!      
 REAL :: rij,rij2
 REAL :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 REAL :: hfacwab,hfacwabi,hfacwabj
 REAL, DIMENSION(ndim) :: dx
!
!  (kernel quantities)
!
 REAL :: q2,q2i,q2j      
 REAL :: wab,wabi,wabj,wabanisoi,wabanisoj,weight
 REAL :: grkern,grkerni,grkernj,grkernanisoi,grkernanisoj
!
!  (grad h terms)
!
 REAL :: dwdhi,dwdhj,dwdhanisoi,dwdhanisoj,dhdrhoi,dhdrhoj
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine density'  
!
!--initialise quantities
!
 nlistdim = ntotal
 ALLOCATE( listneigh(nlistdim) )       ! max size of neighbour list

 DO i=1,npart
    rho(i) = 0.
    gradh(i) = 0.
    gradhaniso(i) = 0.
 ENDDO
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
       hfacwabi = hi1**ndim
       dhdrhoi = -hi*dndim       ! divide by  should use rho(i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
         j = listneigh(n)
         dx(:) = x(:,i) - x(:,j)
         hj = hh(j)
         hj1 = 1./hj
         hj2 = hj*hj
         hfacwabj = hj1**ndim
         dhdrhoj = -hj*dndim       ! see above 
         
         rij2 = DOT_PRODUCT(dx,dx)
         rij = SQRT(rij2)
         q2i = rij2/hi2
         q2j = rij2/hj2       
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!       
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
        IF ((q2i.LT.radkern2).OR.(q2j.LT.radkern2)  &
            .AND. .NOT.(i.GT.npart.AND.j.GT.npart)) THEN
!          PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
!
!--weight self contribution by 1/2
!
           weight = 1.0
           IF (j.EQ.i) weight = 0.5
!       
!--interpolate from kernel table              
!  (use either average h or average kernel gradient)
!
           IF (ikernav.EQ.1) THEN              
!  (using average h)
              hav = 0.5*(hi + hj)
              hav1 = 1./hav
	      hfacwab = hav1**ndim
              q2 = rij2*hav1*hav1
              CALL interpolate_kernel(q2,wab,grkern)
              wab = wab*hfacwab
              grkern = grkern*hfacwab*hav1
	      wabi = wab
	      wabj = wab
           ELSE
              !
              !--calculate both kernels if using the anticlumping term
              !  (so can calculate grad h terms for both kernels)
              !
              if (imhd.ne.0 .and. imagforce.eq.2 .and. ianticlump.eq.1 .and. ikernav.eq.3) then
                 !  (using hi)
                 CALL interpolate_kernels(q2i,wabi,grkerni,wabanisoi,grkernanisoi)
                 wabi = wabi*hfacwabi
                 grkerni = grkerni*hfacwabi*hi1
                 wabanisoi = wabanisoi*hfacwabi
                 grkernanisoi = grkernanisoi*hfacwabi*hi1
                 !  (using hj)
                 CALL interpolate_kernels(q2j,wabj,grkernj,wabanisoj,grkernanisoj)
                 wabj = wabj*hfacwabj
                 grkernj = grkernj*hfacwabj*hj1
                 wabanisoj = wabanisoj*hfacwabj
                 grkernanisoj = grkernanisoj*hfacwabj*hj1
                 !  (calculate average)                
                 !!wab = 0.5*(wabi + wabj)
                 !
                 !--derivative w.r.t. h for grad h correction terms (and dhdrho)
                 !              
                 dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
                 dwdhj = -rij*grkernj*hj1 - ndim*wabj*hj1
                 dwdhanisoi = -rij*grkernanisoi*hi1 - ndim*wabanisoi*hi1
                 dwdhanisoj = -rij*grkernanisoj*hj1 - ndim*wabanisoj*hj1
                 !
                 !--correction term for variable smoothing lengths on anisotropic kernel
                 ! 
                 
                 gradhaniso(i) = gradhaniso(i) + dhdrhoi*pmass(j)*weight*dwdhi
                 gradhaniso(j) = gradhaniso(j) + dhdrhoj*pmass(i)*weight*dwdhj
              else
                 !  (using hi)
                 CALL interpolate_kernel(q2i,wabi,grkerni)
                 wabi = wabi*hfacwabi
                 grkerni = grkerni*hfacwabi*hi1
                 !  (using hj)
                 CALL interpolate_kernel(q2j,wabj,grkernj)
                 wabj = wabj*hfacwabj
                 grkernj = grkernj*hfacwabj*hj1
                 !  (calculate average)                
                 if (ikernav.eq.2) then
		    wab = 0.5*(wabi + wabj)
		    wabi = wab
		    wabj = wab
		 endif
                 !
                 !--derivative w.r.t. h for grad h correction terms (and dhdrho)
                 !              
                 dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
                 dwdhj = -rij*grkernj*hj1 - ndim*wabj*hj1
                 
              endif
                                              
           ENDIF
!
!--calculate density
!
           rho(i) = rho(i) + pmass(j)*wabi*weight
           rho(j) = rho(j) + pmass(i)*wabj*weight
	   
	   IF (ikernav.EQ.3) THEN
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

              gradh(i) = gradh(i) + dhdrhoi*pmass(j)*weight*dwdhi
              gradh(j) = gradh(j) + dhdrhoj*pmass(i)*weight*dwdhj
           ENDIF
!        ELSE
!           PRINT*,' r/h > 2 '      
        
         ENDIF
           
       ENDDO loop_over_neighbours

       iprev = i
       IF (iprev.NE.-1) i = ll(i)              ! possibly should be only IF (iprev.NE.-1)
    ENDDO loop_over_cell_particles
            
 ENDDO loop_over_cells

 RETURN
END SUBROUTINE density
      
