module density_summations
  use get_neighbour_lists
  implicit none
    
contains
!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab (h_a)
!!
!! Also computes the variable smoothing length terms sum_b m_b dW_ab/dh_a
!!
!! This version computes the density on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

  subroutine density(x,pmass,hh,rho,gradh,npart)
    use dimen_mhd
    use debug
    use loguns
 
    use bound
    use kernel
    use linklist
    use options
    use setup_params
    !use matrixcorr
!
!--define local variables
!
    implicit none
    integer, intent(in) :: npart
    real, dimension(:,:), intent(in) :: x
    real, dimension(:), intent(in) :: pmass, hh
    real, dimension(:), intent(out) :: rho, gradh
 
    integer :: i,j,n,idim
    integer :: icell,iprev,nneigh
    integer, dimension(npart) :: listneigh
    integer :: idone
!
!  (particle properties - local copies)
!      
    real :: rij,rij2
    real :: hi,hi1,hav,hav1,hj,hj1,hi2,hj2
    real :: hfacwab,hfacwabi,hfacwabj
    real, dimension(ndim) :: dx
!
!  (kernel quantities)
!
    real :: q2,q2i,q2j      
    real :: wab,wabi,wabj,weight
    real :: grkern,grkerni,grkernj
    real :: dwdhi,dwdhj ! grad h terms
!
!--allow for tracing flow
!      
    if (trace) write(iprint,*) ' Entering subroutine density'  
!
!--initialise quantities
!
    listneigh = 0
    numneigh = 0
    dwdhi = 0.
    dwdhj = 0.

    do i=1,npart
       rho(i) = 0.
       gradh(i) = 0.
       !gradmatrix(:,:,i) = 0.
    enddo
!
!--Loop over all the link-list cells
!
    loop_over_cells: do icell=1,ncellsloop   ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
       call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
       i = ifirstincell(icell)
       idone = -1       ! note density summation includes current particle
       if (i.NE.-1) iprev = i

       loop_over_cell_particles: do while (i.NE.-1)              ! loop over home cell particles

!       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
          idone = idone + 1
          hi = hh(i)
          hi1 = 1./hi
          hi2 = hi*hi
          hfacwabi = hi1**ndim
!
!--for each particle in the current cell, loop over its neighbours
!
          loop_over_neighbours: do n = idone+1,nneigh
             j = listneigh(n)
             dx(:) = x(:,i) - x(:,j)
             hj = hh(j)
             hj1 = 1./hj
             hj2 = hj*hj
             hfacwabj = hj1**ndim
         
             rij2 = dot_product(dx,dx)
             rij = sqrt(rij2)
             q2i = rij2/hi2
             q2j = rij2/hj2       
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!       
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if ((q2i.LT.radkern2).OR.(q2j.LT.radkern2)  &
                  .AND. .not.(i.GT.npart.AND.j.GT.npart)) then
                !if (i.eq.1 .or. j.eq.1) then
                !   print*,' neighbour,r/hi,r/hj,hi,hj:',i,j,sqrt(q2i),sqrt(q2j),hi,hj
                !endif

                if (i.LE.npart) numneigh(i) = numneigh(i) + 1
                if (j.LE.npart .and. j.ne.i) numneigh(j) = numneigh(j) + 1
!
!--weight self contribution by 1/2
!
                weight = 1.0
                if (j.EQ.i) weight = 0.5
!       
!--interpolate from kernel table              
!  (use either average h or average kernel gradient)
!
                if (ikernav.EQ.1) then              
!  (using average h)
                   hav = 0.5*(hi + hj)
                   hav1 = 1./hav
                   hfacwab = hav1**ndim
                   q2 = rij2*hav1*hav1
                   call interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hav1
                   wabi = wab
                   wabj = wab
                   grkerni = grkern
                   grkernj = grkern
                else
              !
              !--calculate both kernels if using the anticlumping term
              !  (so can calculate grad h terms for both kernels)
              !
              !  (using hi)
                   call interpolate_kernel(q2i,wabi,grkerni)
                   wabi = wabi*hfacwabi
                   grkerni = grkerni*hfacwabi*hi1
              !  (using hj)
                   call interpolate_kernel(q2j,wabj,grkernj)
                   wabj = wabj*hfacwabj
                   grkernj = grkernj*hfacwabj*hj1
              !  (calculate average)                
                   if (ikernav.eq.2) then
                      wab = 0.5*(wabi + wabj)
                      wabi = wab
                      wabj = wab
                      grkern = 0.5*(grkerni + grkernj)
                      grkerni = grkern
                      grkernj = grkern
                   endif
              !
              !--derivative w.r.t. h for grad h correction terms (and dhdrho)
              !              
                   dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
                   dwdhj = -rij*grkernj*hj1 - ndim*wabj*hj1
                                              
                endif
!
!--calculate density
!
                rho(i) = rho(i) + pmass(j)*wabi*weight
                rho(j) = rho(j) + pmass(i)*wabj*weight
           
                if (ikernav.EQ.3) then
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

                   gradh(i) = gradh(i) + pmass(j)*weight*dwdhi
                   gradh(j) = gradh(j) + pmass(i)*weight*dwdhj
                endif
                
                !if (i.ne.j) then
                !do idim=1,ndim
                !   gradmatrix(:,idim,i) = gradmatrix(:,idim,i) &
                !               + 2.*pmass(j)*(dx(:))*dx(idim)/rij*grkerni
                !   gradmatrix(:,idim,j) = gradmatrix(:,idim,j) &
                !               + 2.*pmass(i)*(dx(:))*dx(idim)/rij*grkernj
                !enddo
                !endif

!        ELSE
!           PRINT*,' r/h > 2 '      
        
             endif
           
          enddo loop_over_neighbours

          iprev = i
          if (iprev.NE.-1) i = ll(i)  ! possibly should be only IF (iprev.NE.-1)
       enddo loop_over_cell_particles
            
    enddo loop_over_cells
    
    !print*,'end of density, rho, gradh = ',rho(1),gradh(1),hh(1),numneigh(1)
    !print*,'maximum number of neighbours = ',MAXVAL(numneigh),MAXLOC(numneigh),rho(MAXLOC(numneigh))
    !print*,'minimum number of neighbours = ', &
    !       MINVAL(numneigh(1:npart)),MINLOC(numneigh(1:npart)), &
    !       rho(MINLOC(numneigh(1:npart)))

    return
  end subroutine density
      

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
  
  subroutine density_partial(x,pmass,hh,rho,gradh,npart,nlist,ipartlist)
    use dimen_mhd
    use debug
    use loguns
 
    use bound
    use kernel
    use linklist
    !use matrixcorr
!
!--define local variables
!
    implicit none
    integer, intent(in) :: npart
    real, dimension(:,:), intent(in) :: x
    real, dimension(:), intent(in) :: pmass, hh
    real, dimension(:), intent(out) :: rho, gradh
    integer, intent(in) :: nlist
    integer, intent(in), dimension(:) :: ipartlist

    integer :: i,j,n,idim
    integer :: icell,ipart,nneigh,minneigh,minpart
    integer, dimension(npart) :: listneigh
    integer :: icellprev
!
!  (particle properties - local copies)
!      
    real :: rij,rij2
    real :: hi,hi1,hi2
    real :: hfacwabi,hfacgrkerni
    real, dimension(ndim) :: dx
!
!  (kernel quantities)
!
    real :: q2i      
    real :: wabi,grkerni    
    real :: dwdhi ! grad h terms
!
!--allow for tracing flow
!      
    if (trace) write(iprint,*) ' Entering subroutine density_partial'  
!
!--initialise quantities
!
    listneigh = 0

    do ipart=1,nlist
       i = ipartlist(ipart)
       rho(i) = 0.
       gradh(i) = 0.
!       gradmatrix(:,:,i) = 0.
       numneigh(i) = 0
    enddo
    icellprev = 0
!
!--Loop over all the particles in the density list
!
    loop_over_particles: do ipart=1,nlist            ! step through all cells

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
       if (icell.NE.icellprev) then
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       endif
       icellprev = icell

   !!!PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       hi = hh(i)
       hi2 = hi*hi
       hi1 = 1./hi
       
       hfacwabi = hi1**ndim
       hfacgrkerni = hfacwabi*hi1
!
!--loop over current particle's neighbours
!
       loop_over_neighbours: do n = 1,nneigh
          j = listneigh(n)
          dx(:) = x(:,i) - x(:,j)
!
!--calculate averages of smoothing length if using this averaging
!                           
          rij2 = dot_product(dx,dx)
          rij = sqrt(rij2)
          q2i = rij2/hi2
!      
!--do interaction if r/h < compact support size
!
          if (q2i.LT.radkern2) then
             !!!if (i.eq.416) PRINT*,' neighbour,r/h,hi ',j,SQRT(q2i),hi
             numneigh(i) = numneigh(i) + 1
!      
!--interpolate from kernel table (using hi)
!
             call interpolate_kernel(q2i,wabi,grkerni)
             wabi = wabi*hfacwabi
             grkerni = grkerni*hfacgrkerni
!
!--derivative w.r.t. h for grad h correction terms (and dhdrho)
!             
             dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
!
!--calculate density
!
             rho(i) = rho(i) + pmass(j)*wabi
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

             gradh(i) = gradh(i) + pmass(j)*dwdhi
          
             !do idim=1,ndim
             !   gradmatrix(:,idim,i) = gradmatrix(:,idim,i) &
             !               + 2.*pmass(j)*(dx(:))*dx(idim)/rij*grkerni
             !enddo
       
          endif
          
       enddo loop_over_neighbours

    enddo loop_over_particles

    !print*,'finished density_partial, rho, gradh, h =',rho(1),gradh(1),hh(1),numneigh(1)
    !print*,'maximum number of neighbours = ',MAXVAL(numneigh),MAXLOC(numneigh),rho(MAXLOC(numneigh))
    !print*,'minimum number of neighbours = ',MINVAL(numneigh(1:npart)), &
    !       MINLOC(numneigh(1:npart)),rho(MINLOC(numneigh(1:npart)))

    !minneigh = 100000
    !minpart = 1
    !do i=1,nlist
    !   j = ipartlist(i)
    !   if (numneigh(j).lt.minneigh) then
    !      minneigh = numneigh(j)
    !      minpart = j
    !   endif
    !enddo
    !print*,'minimum number of neighbours = ',minneigh,minpart,rho(minpart)

    return
  end subroutine density_partial
      
  
end module density_summations
