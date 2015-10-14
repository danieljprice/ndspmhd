!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------
! Computes an SPH estimate of the second derivatives of a vector quantity
! with various possible operators
!
! Input:
!    iderivtype - type of second derivative operator (see below)
!    npart - number of particles
!    x     - positions
!    pmass - particle masses
!    rho   - density
!    hh    - smoothing length
!    v(:)  - vector quantity to take curl of
!
! Output:
!    del2v(:)    - Laplacian of v
!    graddivv(:) - Gradient of the divergence of v
!
! iderivtype = 1 (default): Brookshaw/Espanol-Revenga style derivatives
! 
!    del2v = \sum m_j/rho_j (v_i - v_j) * 2.0*abs(\nabla W_ij)/rij (h_i)
!    graddivv = \sum m_j/rho_j (5 vij.rij - vij) abs(\nabla W_ij)/rij (h_i)
!
! iderivtype = 2: "standard" second derivatives with kernel
! 
!    del2v = \sum m_j/rho_j (v_j - v_i) nabla^2 W_ij (hi)
!    graddivv = \sum m_j/rho_j ((v_j - v_i).\nabla) nabla W_ij
!
!------------------------------------------------------------------------
module get2ndderivs
 implicit none
 
contains

subroutine get_2ndderivs(iderivtype,npart,x,pmass,rho,hh,v,del2v,graddivv)
 use dimen_mhd,           only:ndim,ndimV,idim
 use debug,               only:trace
 use loguns,              only:iprint
 use kernels,             only:interpolate_kernel_curl,radkern2
 use linklist,            only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use part,                only:ntotal
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart,iderivtype
 real, dimension(ndim,idim),  intent(in) :: x
 real, dimension(idim),       intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(in) :: v
 real, dimension(ndimV,idim), intent(out) :: del2v, graddivv

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,dv,vi,del2vi,graddivvi
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: rho1i,rho1j,del2Wi,del2Wj,rij1,projv
 real, dimension(ntotal) :: h1
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_curl' 
!
!--initialise quantities
!
 listneigh = 0
 dr(:) = 0.
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
 del2v(:,:) = 0.
 graddivv(:,:) = 0.
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

       idone = idone + 1
       hi1 = h1(i)
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
       rho1i  = 1./rho(i)
       vi(:) = v(:,i)
       del2vi(:)    = 0.
       graddivvi(:) = 0.
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself
             dx(:) = x(:,i) - x(:,j)
             hj1 = h1(j)
             
             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
             q2j = rij2*hj1*hj1
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if ((q2i.lt.radkern2).or.(q2j.lt.radkern2) .and. rij2 > 0.) then

                hfacwabj = hj1**ndim
                rij = sqrt(rij2)
                dr(1:ndim) = dx(1:ndim)/(rij + tiny(rij))  ! unit vector
                rij1 = 1./rij
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
!--calculate 2nd derivatives of v
!
                rho1j = 1./rho(j)
                dv = vi(:) - v(:,j)
                projv = dot_product(dv,dr)
                select case(iderivtype)
                case(2)  ! "direct kernel method"
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1
                   grgrkerni = grgrkerni*hfacwabi*hi1*hi1
                   grgrkernj = grgrkernj*hfacwabj*hj1*hj1
                   del2Wi = grgrkerni + (ndim-1)*grkerni*rij1
                   del2Wj = grgrkernj + (ndim-1)*grkernj*rij1

                   del2vi(:)  = del2vi(:)  - pmass(j)*rho1j*dv(:)*del2Wi
                   del2v(:,j) = del2v(:,j) + pmass(i)*rho1i*dv(:)*del2Wj
                   
                   graddivvi(:)  = graddivvi(:)  - pmass(j)*rho1j*(dr(:)*projv*grgrkerni + (dv(:) - projv*dr(:))*grkerni*rij1)
                   graddivv(:,j) = graddivv(:,j) + pmass(i)*rho1i*(dr(:)*projv*grgrkernj + (dv(:) - projv*dr(:))*grkernj*rij1)
                
                case default  ! "Brookshaw method"
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1
                   
                   del2Wi = -2.*grkerni*rij1
                   del2Wj = -2.*grkernj*rij1

                   del2vi(:)  = del2vi(:)  - pmass(j)*rho1j*dv(:)*del2Wi
                   del2v(:,j) = del2v(:,j) + pmass(i)*rho1i*dv(:)*del2Wj
                   
                   graddivvi(:)  = graddivvi(:)  + pmass(j)*rho1j*((ndim + 2)*projv*dr(:) - dv(:))*grkerni*rij1
                   graddivv(:,j) = graddivv(:,j) - pmass(i)*rho1i*((ndim + 2)*projv*dr(:) - dv(:))*grkernj*rij1
                end select

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       del2v(:,i) = del2v(:,i) + del2vi(:)
       graddivv(:,i) = graddivv(:,i) + graddivvi(:)
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 return
end subroutine get_2ndderivs
      
end module get2ndderivs
