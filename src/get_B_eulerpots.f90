!!------------------------------------------------------------------------
!! Computes B = \nabla \alpha_k x \nabla \beta^k
!! (required for the Generalised Euler Potentials)
!! for all particles, doing each pairwise interaction once
!!------------------------------------------------------------------------
module getBeulerpots
 implicit none
 
contains

!
!--iderivtype:
!  1: compute B = \nabla\alpha x \nabla\beta
!  2: compute B = curl (alpha)
!
subroutine get_B_eulerpots(iderivtype,npart,x,pmass,rho,hh,alphapot,betapot,Bfield,remap)
 use dimen_mhd,    only:ndim,ndimV,idim
 use debug,        only:trace
 use loguns,       only:iprint
 
 use kernels,      only:interpolate_kernel_curl,radkern2
 use linklist,     only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms,       only:gradh
 use setup_params, only:hfact
 use part,         only:itype,ntotal
 use matrixcorr,   only:dxdx,ndxdx,idxdx,jdxdx
 use bound,        only:ireal
 use options,      only:iuse_exact_derivs
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart,iderivtype
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(inout) :: alphapot
 real, dimension(ndimV,idim), intent(inout)  :: betapot
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
 real, dimension(ndimV,ndimV) :: gradalphai,gradbetai
 real, dimension(ndxdx) :: dxdxi
 real, dimension(6) :: rmatrix
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: rho21i,rho21gradhi,denom,ddenom
 real, dimension(ntotal) :: h1
 real, dimension(ndimV,ndimV,ntotal) :: gradalpha, gradbeta
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_B_eulerpots' 
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
 dxdx(:,:) = 0.
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
       dxdxi(:) = 0.
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

                   do k=1,ndxdx
                      dxdxi(k) = dxdxi(k) - pmass(j)*(dx(idxdx(k)))*dr(jdxdx(k))*grkerni
                      dxdx(k,j) = dxdx(k,j) - pmass(i)*(dx(idxdx(k)))*dr(jdxdx(k))*grkernj
                   enddo

                end select

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       gradalpha(:,:,i) = gradalpha(:,:,i) + gradalphai(:,:)
       gradbeta(:,:,i) = gradbeta(:,:,i) + gradbetai(:,:)
       dxdx(:,i) = dxdx(:,i) + dxdxi(:)
       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    select case(iderivtype)
    case default
       if (iuse_exact_derivs.gt.0) then
          call compute_rmatrix(dxdx(:,i),rmatrix,denom,ndim)
          !print*,i,' normal derivs, gradalpha^z = ',gradalpha(:,3,i)*gradh(i)/rho(i)
          if (abs(denom).gt.epsilon(denom)) then
             ddenom = 1./denom
             do k=1,ndimV
                call exactlinear(gradalpha(:,k,i),gradalpha(:,k,i),rmatrix,ddenom)
                call exactlinear(gradbeta(:,k,i),gradbeta(:,k,i),rmatrix,ddenom)
             enddo
             !print*,i,' exact derivs, gradalpha^z = ',gradalpha(:,3,i)
          else
             print*,'WARNING: denominator collapsed in exact linear deriv = ',denom
             gradalpha(:,:,i) = gradalpha(:,:,i)*gradh(i)/rho(i)
             gradbeta(:,:,i) = gradbeta(:,:,i)*gradh(i)/rho(i)
          endif
          !read*
       else
          gradalpha(:,:,i) = gradalpha(:,:,i)*gradh(i)/rho(i)
          gradbeta(:,:,i) = gradbeta(:,:,i)*gradh(i)/rho(i)
       endif

       if (iderivtype.eq.2) then
          !--compute B = curl alpha (i.e., eps_ijk grad^j alpha^k)
          call curl3D_epsijk(gradalpha(:,:,i),Bfield(:,i))
          !print*,i,'Bfield from curl A= ',Bk(:)
       else
          !--compute B = grad alpha x grad beta
          do k=1,ndimV
             call cross_product3D(gradalpha(:,k,i),gradbeta(:,k,i),Bk)
             !print*,i,' B = B + ',Bk(:)
             Bfield(:,i) = Bfield(:,i) + Bk(:)
          enddo
          !print*,i,' Bfield from grad alpha x grad beta = ',Bfield(:,i)
          !read*
          !--remap to a new set of potentials
          if (remap) then
             !print*,i,' old vector potential = ',alphapot(:,i)
             alphapoti(:) = 0.
             do k=1,ndimV
                alphapoti(:) = alphapoti(:) + alphapot(k,i)*gradbeta(:,k,i)
             enddo
             !print*,i,' new vector potential = ',alphapoti(:)
             alphapot(:,i) = alphapoti(:)
             betapot(1:ndim,i) = x(1:ndim,i)
          endif
       endif
    end select
 enddo

 return
end subroutine get_B_eulerpots

!----------------------------------------------------------------
!+
!  Computes matrix terms needed for exact linear derivatives
!+
!----------------------------------------------------------------
subroutine compute_rmatrix(dxdxi,rmatrix,denom,ndim)
 !use matrixcorr, only:ndxdx
 implicit none
 real, dimension(:), intent(in) :: dxdxi
 real, dimension(6), intent(out) :: rmatrix
 real, intent(out) :: denom
 integer, intent(in) :: ndim
 real :: rxxi,rxyi,rxzi,ryyi,ryzi,rzzi

 !print*,'ndim = ',ndim,'got dxdx = ',dxdxi,' ndxdx=',ndxdx
 
 rxxi = dxdxi(1)
 if (ndim.ge.2) then
    rxyi = dxdxi(2)
    ryyi = dxdxi(3)
 else
    rxyi = 0.
    ryyi = 1.
 endif
 if (ndim.ge.3) then
    rxzi = dxdxi(4)
    ryzi = dxdxi(5)
    rzzi = dxdxi(6)
 else
    rxzi = 0.
    ryzi = 0.
    rzzi = 1.
 endif
 
 !print*,' got dxdx = ',dxdxi(:)

 denom = rxxi*ryyi*rzzi + 2.*rxyi*rxzi*ryzi &
       - rxxi*ryzi*ryzi - ryyi*rxzi*rxzi - rzzi*rxyi*rxyi

 !print*,' gives denom = ',rxxi*ryyi*rzzi,2.*rxyi*rxzi*ryzi,-rxxi*ryzi*ryzi,-ryyi*rxzi*rxzi,-rzzi*rxyi*rxyi

 rmatrix(1) = ryyi*rzzi - ryzi*ryzi    ! xx
 rmatrix(2) = rxzi*ryzi - rzzi*rxyi    ! xy
 rmatrix(3) = rxyi*ryzi - rxzi*ryyi    ! xz
 rmatrix(4) = rzzi*rxxi - rxzi*rxzi    ! yy
 rmatrix(5) = rxyi*rxzi - rxxi*ryzi    ! yz
 rmatrix(6) = rxxi*ryyi - rxyi*rxyi    ! zz
 !print*,' rmatrix = ',rmatrix(:), 'denom =',denom
 !read*
end subroutine compute_rmatrix

!----------------------------------------------------------------
!+
!  Internal subroutine that inverts the matrix to get an
!  exact linear derivative
!+
!----------------------------------------------------------------
pure subroutine exactlinear(gradA,dAin,rmatrix,ddenom)
 implicit none
 real, dimension(3), intent(out) :: gradA
 real, dimension(3), intent(in)  :: dAin
 real, intent(in), dimension(6) :: rmatrix
 real, intent(in)  :: ddenom
 real, dimension(size(dAin)) :: dA
 
 !--make local copy of dA to allow input array (dAin)
 !  to be same as output (gradA)
 dA = dAin

 gradA(1) =(dA(1)*rmatrix(1) + dA(2)*rmatrix(2) + dA(3)*rmatrix(3))*ddenom
 gradA(2) =(dA(1)*rmatrix(2) + dA(2)*rmatrix(4) + dA(3)*rmatrix(5))*ddenom
 gradA(3) =(dA(1)*rmatrix(3) + dA(2)*rmatrix(5) + dA(3)*rmatrix(6))*ddenom

 !--the above is equivalent to:
 !gradAx =(dAx*termxx + dAy*termxy + dAz*termxz)*ddenom
 !gradAy =(dAx*termxy + dAy*termyy + dAz*termyz)*ddenom
 !gradAz =(dAx*termxz + dAy*termyz + dAz*termzz)*ddenom

 return
end subroutine exactlinear

end module getBeulerpots
