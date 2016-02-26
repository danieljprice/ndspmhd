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
! Computes anisotropic diffusion terms in SPH
!
!------------------------------------------------------------------------
module aniso_diffusion
 implicit none
 
contains

subroutine get_aniso_diffusion(npart,x,pmass,rho,hh,u,dudt,k_tensor)
 use dimen_mhd,           only:ndim,ndimV
 use debug,               only:trace
 use loguns,              only:iprint
 use kernels,             only:interpolate_kernel,radkern2
 use linklist,            only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use part,                only:ntotal
 use utils,               only:delta_fn
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart
 real, dimension(:,:),  intent(in) :: x
 real, dimension(:),    intent(in) :: pmass,rho,hh, u
 real, dimension(:),    intent(out) :: dudt
 real, dimension(ndim,ndim), intent(in) :: k_tensor

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr
 real :: q2i,q2j,wi,wj,grkerni,grkernj
 real :: rho1i,rho1j,rij1,rhoi,rhoj,term,grkern
 real, dimension(ntotal) :: h1
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_aniso_diffusion' 
!
!--initialise quantities
!
 listneigh = 0
 dr(:) = 0.
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
 dudt(:) = 0.
 print*,' k_tensor = ',k_tensor
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
       rhoi   = rho(i)
       rho1i  = 1./rhoi
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
                call interpolate_kernel(q2i,wi,grkerni)
                !  (using hj)
                call interpolate_kernel(q2j,wj,grkernj)

                rhoj  = rho(j)
                rho1j = 1./rhoj

                grkerni = grkerni*hfacwabi*hi1
                grkernj = grkernj*hfacwabj*hj1
                grkern = 0.5*(grkerni + grkernj)
!
!--general anisotropic diffusion tensor
!
                term = get_term(ndim,k_tensor,dr)
                dudt(i) = dudt(i) + pmass(j)*rho1j*(u(i) - u(j))*term*grkern*rij1
                dudt(j) = dudt(j) - pmassi*rho1i*(u(i) - u(j))*term*grkern*rij1

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells
! print*,' HERE, dtcourant = ',dtcourant
 
 return
end subroutine get_aniso_diffusion

pure real function get_term(ndim,k_tensor,dr)
 use utils, only:delta_fn
 integer, intent(in) :: ndim
 real, dimension(ndim,ndim), intent(in) :: k_tensor
 real, dimension(ndim), intent(in) :: dr
 integer :: i,j
 
 get_term = 0.
 do i=1,ndim
    do j=1,ndim
       get_term = get_term + k_tensor(i,j)*(5.*dr(i)*dr(j) - delta_fn(i,j))
    enddo
 enddo

end function get_term

end module aniso_diffusion
