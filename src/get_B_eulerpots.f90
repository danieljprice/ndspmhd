!!------------------------------------------------------------------------
!! Computes B = \nabla \alpha_k x \nabla \beta^k
!! (required for the Generalised Euler Potentials)
!! for all particles, doing each pairwise interaction once
!!------------------------------------------------------------------------
module getBeulerpots
 implicit none
 
contains

subroutine get_B_eulerpots(iderivtype,npart,x,pmass,rho,hh,alphapot,betapot,Bfield,remap)
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
 integer, intent(in) :: npart,iderivtype
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(inout) :: alphapot
 real, dimension(ndim,idim), intent(inout)  :: betapot
 real, dimension(ndimV,idim), intent(out) :: Bfield
 logical, intent(in) :: remap
 real :: weight

 integer :: i,j,n,k
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,dalpha,dbeta,Bk,alphapoti,betapoti
 real, dimension(ndim,ndimV) :: gradalphai,gradbetai
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: rho21i,rho21gradhi
 real, dimension(ntotal) :: h1
 real, dimension(ndim,ndimV,ntotal) :: gradalpha, gradbeta
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_curl' 
!
!--initialise quantities
!
 listneigh = 0
 gradalpha = 0.
 gradbeta = 0.
 dr(:) = 0.
 weight = 1./hfact**ndim
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
 Bfield = 0.
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
       alphapoti(:) = alphapot(:,i)
       betapoti(:) = betapot(:,i)
       rho21i = 1./rho(i)**2
       rho21gradhi = rho21i*gradh(i)
       gradalphai(:,:) = 0.
       gradbetai(:,:) = 0.
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
                dr(1:ndim) = dx(1:ndim)/(rij + tiny(rij))  ! unit vector
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                !  (using hi)
                call interpolate_kernel_curl(q2i,grkerni,grgrkerni)
                !  (using hj)
                call interpolate_kernel_curl(q2j,grkernj,grgrkernj)
!
!--calculate grad alpha_k and grad beta^k
!
                select case(iderivtype)
                case default  ! default mass weighted grad
                              ! i.e. grad A = 1/(omega rhoi) \sum mj (A_j - A_i) grad W
                              ! -- divide by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   dalpha(1:ndimV) = alphapot(1:ndimV,j) - alphapoti(1:ndimV)
                   dbeta(1:ndimV) = betapot(1:ndimV,j) - betapoti(1:ndimV)
!                   call cross_product3D(dB,dr,curlBterm)

                   do k = 1,ndimV
                      gradalphai(:,k) = gradalphai(:,k) + pmass(j)*dalpha(k)*dr(:)*grkerni
                      gradalpha(:,k,j) = gradalpha(:,k,j) + pmassi*dalpha(k)*dr(:)*grkernj
                      gradbetai(:,k) = gradbetai(:,k) + pmass(j)*dbeta(k)*dr(:)*grkerni
                      gradbeta(:,k,j) = gradbeta(:,k,j) + pmassi*dbeta(k)*dr(:)*grkernj
                   enddo

                end select

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       gradalpha(:,:,i) = gradalpha(:,:,i) + gradalphai(:,:)
       gradbeta(:,:,i) = gradbeta(:,:,i) + gradbetai(:,:)
       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    select case(iderivtype)
    case default
       gradalpha(:,:,i) = gradalpha(:,:,i)*gradh(i)/rho(i)
       gradbeta(:,:,i) = gradbeta(:,:,i)*gradh(i)/rho(i)
       
       !--compute B = grad alpha x grad beta
       do k=1,ndimV
          call cross_product3D(gradalpha(:,k,i),gradbeta(:,k,i),Bk)
          !print*,i,' B = B + ',Bk(:)
          Bfield(:,i) = Bfield(:,i) + Bk(:)
       enddo
       !--remap to a new set of potentials
       if (remap) then
          !print*,i,' old vector potential = ',alphapot(:,i)
          alphapoti(:) = 0.
          do k=1,ndimV
             alphapoti(:) = alphapoti(:) + alphapot(k,i)*gradbeta(:,k,i)
          enddo
          !print*,i,' new vector potential = ',alphapoti(:)
          alphapot(:,i) = alphapoti(:)
          betapot(:,i) = x(:,i)
       endif
    end select
 enddo

 return
end subroutine get_B_eulerpots
      
end module getBeulerpots
