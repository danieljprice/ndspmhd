!!------------------------------------------------------------------------
!! Computes dB/dt = eta nabla^2 B implicitly
!!------------------------------------------------------------------------
module resistivity
 implicit none
 
contains

!
!--does implicit B diffusion
!
subroutine Bdiffusion(npart,x,pmass,rho,hh,Bfield,dBevoldt,dt)
 use dimen_mhd,    only:ndim,ndimV,idim
 use debug,        only:trace
 use loguns,       only:iprint
 
 use kernels,      only:interpolate_kernel,radkern2
 use linklist,     only:ll,ifirstincell,ncellsloop,numneigh
 use get_neighbour_lists, only:get_neighbour_list_partial
 use hterms,       only:gradh
 use setup_params, only:hfact
 use part,         only:itype,ntotal,rho0
 use bound,        only:ireal
 use options,      only:ibound,etamhd
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(in) :: Bfield
 real, dimension(ndimV,idim), intent(inout) :: dBevoldt
 real, intent(in) :: dt

 integer :: i,j,n,isweep,nsweeps
 integer :: icell,iprev,nneigh,nneighi
 integer, dimension(npart) :: listneigh ! neighbour list
 real :: rij,rij2
 real :: hi1,hi21
 real :: hfacwabi,rho1j,sumdenom
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,Bi,dB,sumB,sumBin
 real :: q2i,grkerni,wabi,dB2,dBmax,term
 real, dimension(ndimV,ntotal) :: Bfieldnew
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine Bdiffusion' 
!
!--initialise quantities
!
 listneigh = 0
 dr(:) = 0.
 do i=1,ntotal
    Bfieldnew(:,i) = Bfield(:,i)
 enddo
 nsweeps = 10
 
 sweeps: do isweep = 1,nsweeps 

!
!--update Bfieldnew on ghosts
!
 if (any(ibound.gt.1)) then
    do i=npart+1,ntotal
       j = ireal(i)
       Bfieldnew(:,i) = Bfieldnew(:,j)
    enddo
 endif

 dBmax = 0.
!
!--loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list_partial(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1)          ! loop over home cell particles

       !       print*,'doing particle ',i,nneigh,' neighbours',hh(i)
       !!hi = hh(i)
       if (i.gt.npart) then
          i = ll(i)
          cycle loop_over_cell_particles
       endif
       hi1 = 1./hh(i)
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       nneighi = 0
       sumdenom = 0.
       sumB(:) = 0.
       sumBin(:) = 0.
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = 1,nneigh
          j = listneigh(n)
          if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself
             dx(:) = x(:,i) - x(:,j)
             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if (q2i.lt.radkern2) then
                rij = sqrt(rij2)
                dr(1:ndim) = dx(1:ndim)/(rij + tiny(rij))  ! unit vector
                rho1j = 1./rho(j)
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                !  (using hi)
                call interpolate_kernel(q2i,wabi,grkerni)
!
!--calculate nabla^2 B
!
                grkerni = grkerni*hfacwabi*hi1
                
                term = -2.*etamhd*pmass(j)*rho1j*grkerni/rij
                sumdenom = sumdenom + term
                sumB(:) = sumB(:) + term*Bfieldnew(:,j)
                sumBin(:) = sumBin(:) + term*Bfield(:,j)
                nneighi = nneighi + 1

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours

       Bi(:) = Bfieldnew(:,i)
 
       if (nsweeps.eq.1) then
          !--explicit
          Bfieldnew(:,i) = Bfield(:,i) - dt*(Bfield(:,i)*sumdenom - sumBin(:))
       else
          !--implicit
          Bfieldnew(:,i) = (Bfield(:,i) + dt*sumB(:))/(1. + dt*sumdenom)
       endif
       !print*,' particle ',i,' numneigh = ',numneigh(i),nneighi
       
       dB = Bfieldnew(:,i) - Bi(:)
       dB2 = dot_product(dB,dB)
       dBmax = max(dBmax,dB2)

       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 !print*,' end of sweep ',isweep,' dBmax = ',dBmax,' dt =  ',dt

 enddo sweeps

 if (dt.gt.0.) then
 do i=1,npart
    dBevoldt(:,i) = dBevoldt(:,i) + (Bfieldnew(:,i) - Bfield(:,i))/dt
 enddo
 endif
 !read*

 return
end subroutine Bdiffusion

end module resistivity
