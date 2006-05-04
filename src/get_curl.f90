!!------------------------------------------------------------------------
!! Computes an SPH estimate of the curl of a vector quantity
!! This version computes the curl for all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------
module getcurl
 implicit none
 
contains

subroutine get_curl(npart,x,pmass,rho,hh,Bvec,curlB)
 use dimen_mhd, only:ndim,ndimV,idim
 use debug, only:trace
 use loguns, only:iprint
 
 use kernels, only:interpolate_kernel,radkern2
 use linklist, only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms, only:gradh
 use setup_params, only:hfact
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(in) :: Bvec
 real, dimension(ndimV,idim), intent(out) :: curlB
 real :: weight

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi,hi1,hj,hj1,hi2,hj2
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr, dB, Bi, curlBi
 real :: q2i,q2j,wabi,wabj,grkerni,grkernj
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_curl' 
!
!--initialise quantities
!
 listneigh = 0
 curlB = 0.
 dr(:) = 0.
 weight = 1./hfact**ndim
!
!--loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1)          ! loop over home cell particles

       !       print*,'doing particle ',i,nneigh,' neighbours',hh(i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
       Bi(:) = Bvec(:,i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             hj = hh(j)
             hj1 = 1./hj
             hj2 = hj*hj
             
             rij2 = dot_product(dx,dx)
             q2i = rij2/hi2
             q2j = rij2/hj2     
             !          print*,' neighbour,r/h,dx,hi,hj ',j,sqrt(q2),dx,hi,hj
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then

                hfacwabj = hj1**ndim
                rij = sqrt(rij2)
                dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                !  (using hi)
                call interpolate_kernel(q2i,wabi,grkerni)
!                wabi = wabi*hfacwabi
                grkerni = grkerni*hi1 !!*hfacwabi*hi1
                !  (using hj)
                call interpolate_kernel(q2j,wabj,grkernj)
!                wabj = wabj*hfacwabj
                grkernj = grkernj*hj1 !!*hfacwabj*hj1
!
!--calculate curl of Bvec (NB dB is 3-dimensional, dr is also but zero in parts)
!
                dB(1:ndimV) = Bi(1:ndimV) - Bvec(1:ndimV,j)
!                if (ndim.eq.3) then
                   call cross_product3D(dB,dr,curlBi)
!                   curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
!                   curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
!                   curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
!                elseif (ndim.eq.2) then  ! just Az in 2D
!                   curlBi = 0.
!                   curlBi(1) = -dB(1)*dr(2) ! replace dB(3) by dB(1)
!                   curlBi(2) = dB(1)*dr(1)
!                endif
                !
                !--compute rho * current density j
                !
                curlB(:,i) = curlB(:,i) + pmass(j)*curlBi(:)*grkerni*hfacwabi
                curlB(:,j) = curlB(:,j) + pmassi*curlBi(:)*grkernj*hfacwabj
!                curlB(:,i) = curlB(:,i) + curlBi(:)*grkerni
!                curlB(:,j) = curlB(:,j) + curlBi(:)*grkernj
                !!print*,'weight = ',weight,' m/rho h^3 = ',pmass(j)/rho(j)*hfacwabj
                
             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    curlB(:,i) = curlB(:,i)*gradh(i)/rho(i)
!    curlB(:,i) = weight*curlB(:,i)*gradh(i)
!    print*,i,curlB(:,i)
 enddo

 return
end subroutine get_curl
      
end module getcurl
