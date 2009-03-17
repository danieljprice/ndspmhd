!!------------------------------------------------------------------------
!! Computes an SPH estimate of the curl of a vector quantity
!! This version computes the curl for all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------
module getcurl
 implicit none
 
contains

subroutine get_curl(npart,x,pmass,rho,hh,Bvec,curlB,curlBgradh)
 use dimen_mhd, only:ndim,ndimV,idim
 use debug, only:trace
 use loguns, only:iprint
 
 use kernels, only:interpolate_kernel_curl,radkern2
 use linklist, only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms, only:gradh
 use setup_params, only:hfact
 use part, only:itype,ntotal
 use bound, only:ireal
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(in) :: Bvec
 real, dimension(ndimV,idim), intent(out) :: curlB
 real, dimension(ndimV,idim), intent(out), optional :: curlBgradh
 real :: weight

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi,hi1,hj,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr, dB, Bi, curlBi, curlBterm, curlBgradhi
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: dgradwdhi,dgradwdhj
 real, dimension(ntotal) :: h1
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
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
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
       !!hi = hh(i)
       hi1 = h1(i) !!1./hi
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
       Bi(:) = Bvec(:,i)
       curlBi(:) = 0.
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             !!hj = hh(j)
             hj1 = h1(j) !!1./hj
             
             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
             q2j = rij2*hj1*hj1    
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
                call interpolate_kernel_curl(q2i,grkerni,grgrkerni)
!                wabi = wabi*hfacwabi
                grkerni = grkerni*hi1 !!*hfacwabi*hi1
                dgradwdhi = -(ndim+1.)*grkerni*hfacwabi*hi21 - rij*hi21*grgrkerni*hfacwabi
                !  (using hj)
                call interpolate_kernel_curl(q2j,grkernj,grgrkernj)
!                wabj = wabj*hfacwabj
                grkernj = grkernj*hj1 !!*hfacwabj*hj1
                dgradwdhj = -(ndim+1.)*grkernj*hfacwabj*hj1*hj1 - rij*hj1**2*grgrkernj*hfacwabj
!
!--calculate curl of Bvec (NB dB is 3-dimensional, dr is also but zero in parts)
!
!--form 1 (dB using weights) -- multiply by weight below
                dB(1:ndimV) = Bi(1:ndimV) - Bvec(1:ndimV,j)
                call cross_product3D(dB,dr,curlBterm)
!                curlBi(:) = curlBi(:) + curlBterm(:)*grkerni
!                curlB(:,j) = curlB(:,j) + curlBterm(:)*grkernj                

!--form 2 (dB, m_j/rho_i with gradh) -- divide by rho(i) below
                curlBi(:) = curlBi(:) + pmass(j)*curlBterm(:)*grkerni*hfacwabi
                curlB(:,j) = curlB(:,j) + pmassi*curlBterm(:)*grkernj*hfacwabj

                if (present(curlBgradh)) then
                curlBgradhi(:) = curlBgradhi(:) + pmass(j)*curlBterm(:)*dgradwdhi
                curlBgradh(:,j) = curlBgradh(:,j) + pmassi*curlBterm(:)*dgradwdhj
                endif

!--form 3 (m_j/rho_j) either with dB (use curl above) or with just B (uncomment curls below)
!                call cross_product3D(-Bvec(1:ndimV,j),dr,curlBterm)
!                curlBi(:) = curlBi(:) + pmass(j)/rho(j)*curlBterm(:)*grkernj*hfacwabj
!                call cross_product3D(Bi(1:ndimV),dr,curlBterm)
!                curlB(:,j) = curlB(:,j) + pmassi/rho(i)*curlBterm(:)*grkerni*hfacwabi

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       curlB(:,i) = curlB(:,i) + curlBi(:)
       curlBgradh(:,i) = curlBgradh(:,i) + curlBgradhi(:)
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    curlB(:,i) = curlB(:,i)*gradh(i)/rho(i)
!    curlB(:,i) = weight*curlB(:,i) !!*gradh(i)
!    curlB(:,i) = rho(i)*curlB(:,i) !!*gradh(i)
 enddo

 return
end subroutine get_curl
      
end module getcurl
