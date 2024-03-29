!!------------------------------------------------------------------------
!! Computes B = \nabla \alpha_k x \nabla \x0^k
!! (required for the Generalised Euler Potentials)
!! for all particles, doing each pairwise interaction once
!!------------------------------------------------------------------------
module getBeulerpots
 implicit none
 
contains

!
!--iderivtype:
!  1: compute B = \nabla\alpha x \nabla\x0
!  2: compute B = curl (alpha)
!  3: compute B.grad x0 = Bevol (i.e, B/rho.grad X^0 = B/rho_0)
!
subroutine get_B_eulerpots(iderivtype,npart,x,pmass,rho,hh,alphapot,x0,Bfield,remap,remap_tol)
 use dimen_mhd,    only:ndim,ndimV,idim
 use debug,        only:trace
 use loguns,       only:iprint
 
 use kernels,      only:interpolate_kernel_curl,radkern2
 use linklist,     only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms,       only:gradh
 use setup_params, only:hfact
 use part,         only:itype,ntotal,rho0
 use matrixcorr,   only:dxdx,ndxdx,idxdx,jdxdx
 use bound,        only:ireal
 use options,      only:iuse_exact_derivs,ibound
 use utils,        only:cross_product3D, curl3D_epsijk, matrixinvert3D, det
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart,iderivtype
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(inout) :: alphapot
 real, dimension(ndimV,idim), intent(inout)  :: x0
 real, dimension(ndimV,idim), intent(out) :: Bfield
 logical, intent(inout) :: remap
 real, intent(in), optional :: remap_tol
 real :: weight

 integer :: i,j,n,k,ierr
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone,ierrmax,iJmax
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi,detJ1,detJ2,err,errtot,errmax,detJmax
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,dalpha,dx0,Bk,alphapoti,x0i
 real, dimension(ndimV,ndimV) :: gradalphai,gradx0i
 real, dimension(ndxdx) :: dxdxi
 real, dimension(6) :: rmatrix
 real, dimension(3,3) :: dxdx0i
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: rho21i,rho21gradhi,denom,ddenom,detJmin
 real, dimension(ntotal) :: h1
 real, dimension(ndimV,ndimV,ntotal) :: gradalpha, gradx0
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_B_eulerpots' 
!
!--initialise quantities
!
 listneigh = 0
 gradalpha = 0.
 gradx0 = 0.
 dr(:) = 0.
 weight = 1./hfact**ndim
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
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
       x0i(:) = x0(:,i)
       rho21i = 1./rho(i)**2
       rho21gradhi = rho21i*gradh(i)
       gradalphai(:,:) = 0.
       gradx0i(:,:) = 0.
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
!--calculate grad alpha_k and grad x0^k
!
                select case(iderivtype)
                case default  ! default mass weighted grad
                              ! i.e. grad A = 1/(omega rhoi) \sum mj (A_j - A_i) grad W
                              ! -- divide by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   dalpha(1:ndimV) = alphapot(1:ndimV,j) - alphapoti(1:ndimV)
                   dx0(1:ndimV) = x0(1:ndimV,j) - x0i(1:ndimV)
!                   call cross_product3D(dB,dr,curlBterm)

                   do k = 1,ndimV
                      gradalphai(:,k) = gradalphai(:,k) + pmass(j)*dalpha(k)*dr(:)*grkerni
                      gradalpha(:,k,j) = gradalpha(:,k,j) + pmassi*dalpha(k)*dr(:)*grkernj
                      gradx0i(:,k) = gradx0i(:,k) + pmass(j)*dx0(k)*dr(:)*grkerni
                      gradx0(:,k,j) = gradx0(:,k,j) + pmassi*dx0(k)*dr(:)*grkernj
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
       gradx0(:,:,i) = gradx0(:,:,i) + gradx0i(:,:)
       dxdx(:,i) = dxdx(:,i) + dxdxi(:)
       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

!
!--finalise the gradient terms
!  (we need to do this in a separate loop because we use the determinant
!   as a criterion for whether or not to remap the potentials)
!
 detJmin = huge(detJmin)
 do i=1,npart
    if (iuse_exact_derivs.gt.0) then
       call compute_rmatrix(dxdx(:,i),rmatrix,denom,ndim)
       !print*,i,' normal derivs, gradalpha^z = ',gradalpha(:,3,i)*gradh(i)/rho(i)
       if (abs(denom).gt.epsilon(denom)) then
          if (denom.lt.0.) print*,'WARNING: determinant < 0 in exact derivs particle ',i
          ddenom = 1./denom
          do k=1,ndimV
             call exactlinear(gradalpha(:,k,i),gradalpha(:,k,i),rmatrix,ddenom)
             call exactlinear(gradx0(:,k,i),gradx0(:,k,i),rmatrix,ddenom)
          enddo
          !print*,i,' exact derivs, gradalpha^z = ',gradalpha(:,3,i)
       else
          print*,'WARNING: denominator collapsed in exact linear deriv = ',denom
          gradalpha(:,:,i) = gradalpha(:,:,i)*gradh(i)/rho(i)
          gradx0(:,:,i) = gradx0(:,:,i)*gradh(i)/rho(i)
       endif
       !read*
    else
       gradalpha(:,:,i) = gradalpha(:,:,i)*gradh(i)/rho(i)
       gradx0(:,:,i) = gradx0(:,:,i)*gradh(i)/rho(i)
    endif
    detJ1 = det(gradx0(:,:,i))
    
    detJmin = min(detJ1,detJmin)
 enddo

 if (present(remap_tol)) then
    remap = .false.
    if (detJmin < remap_tol) then
       print*,' minimum det(J) = ',detJmin, ': REMAPPING...'
       remap = .true.
    else
       print*,' minimum det(J) = ',detJmin
    endif
 endif

 !--determine whether or not to remap the particles
 errtot = 0.
 errmax = 0.
 ierrmax = 0
 detJmax = 0.
 iJmax = 0
 do i=1,npart
    !
    !--determinant of 3x3 matrix dx0/dx
    !
    detJ1 = det(gradx0(:,:,i))
    select case(iderivtype)
    case(3,4)
       !--remapping of B/rho
       !  (alphapot is B/rho_0, send out B/rho new as Bfield)
       alphapoti(:) = alphapot(:,i) ! copy so that Bfield can be same as alphapot (input)
       call matrixinvert3D(gradx0(:,:,i),dxdx0i,ierr)
       !print*,'matrix inverse = ',matmul(gradx0(:,:,i),dxdx0i)
       if (ierr.ne.0) then ! NB, this is non-fatal (just means B = 0)
          print*,'WARNING: could not invert matrix in remapping, det = ',detJ1
       endif
!          if (detJ1.le.0.) then
!             print*,'ERROR: determinant collapsed in remapping, particle ',i
!             print*,'J = ',detJ1, 'should be rho0/rho = ',rho0(i)/rho(i)
!          endif                
       detJ2 = rho(i)/rho0(i)
       err = abs(detJ1-detJ2)
       errtot = errtot + err
       if (err.gt.errmax) then             
          errmax = err
          ierrmax = i
       endif
       if (detJ1.gt.detJmax) then
          detJmax = detJ1
          iJmax = i
       endif
          !print*,' abs error in J vs rho0/rho = ',abs(detJ1-detJ2)
!          print*,'J = det(dx/dx0) = ',detJ1,' rho_0/rho = ',detJ2

       if (iderivtype.eq.4) then
          Bfield(:,i) = matmul(alphapoti,dxdx0i)*detJ1  ! this is B (new)
       else
          Bfield(:,i) = matmul(alphapoti,dxdx0i)  ! this is B/rho (new)
       endif
!          do k=1,ndimV
!             Bfield(k,i) = dot_product(alphapoti(:),dxdx0i(:,k))
!          enddo
       if (remap) then
          alphapot(:,i) = Bfield(:,i)
          x0(1:ndim,i) = x(1:ndim,i)
          rho0(i) = rho(i)
       endif
    case(2)
       !--compute B = curl alpha (i.e., eps_ijk grad^j alpha^k)
       call curl3D_epsijk(gradalpha(:,:,i),Bfield(:,i))
       !print*,i,'Bfield from curl A= ',Bk(:)
    case default
       Bfield(:,i) = 0.
       !--compute B = grad alpha x grad x0
       do k=1,ndimV
          call cross_product3D(gradalpha(:,k,i),gradx0(:,k,i),Bk)
          !print*,i,' B = B + ',Bk(:)
          Bfield(:,i) = Bfield(:,i) + Bk(:)
       enddo
       !print*,i,' Bfield from grad alpha x grad x0 = ',Bfield(:,i)
       !read*
       !--remap to a new set of potentials
       if (remap) then
          !print*,i,' old vector potential = ',alphapot(:,i)
          alphapoti(:) = 0.
          do k=1,ndimV
             alphapoti(:) = alphapoti(:) + alphapot(k,i)*gradx0(:,k,i)
          enddo
          !print*,i,' new vector potential = ',alphapoti(:)
          alphapot(:,i) = alphapoti(:)
          x0(1:ndim,i) = x(1:ndim,i)
       endif
    end select
 enddo

 if (iderivtype.eq.3 .or. iderivtype.eq.4) then
    print*,' remapping errors, max = ',errmax,' mean = ',errtot/real(npart)
    print*,' worst on particle ',ierrmax,'   J = ',det(gradx0(:,:,ierrmax)),&
           'rho/rho0 = ',rho(ierrmax)/rho0(ierrmax)
    call matrixinvert3D(gradx0(:,:,ierrmax),dxdx0i,ierr)
    print*,' inverse on part   ',ierrmax,' 1/J = ',det(dxdx0i),&
           'rho0/rho = ',rho0(ierrmax)/rho(ierrmax)

    print*,' max J on particle ',iJmax,' J = ',det(gradx0(:,:,iJmax)),&
           'rho/rho0 = ',rho(iJmax)/rho0(iJmax)
 endif
 
 if (remap .and. iderivtype.ne.2) then
    !--remap x0 for all particles
    x0(1:ndim,:) = x(1:ndim,:)
    !--copy remapped Bevol onto boundary/ghost particles 
    if (any(ibound.eq.1)) then
       do i=1,npart
          if (itype(i).eq.1) then ! fixed particles
             j = ireal(i)
             alphapot(:,i) = alphapot(:,j)
          endif
       enddo
    endif
    if (any(ibound.gt.1)) then  ! ghost particles
       do i=npart+1,ntotal
          j = ireal(i)
          alphapot(:,i) = alphapot(:,j)
       enddo
    endif
 endif

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
