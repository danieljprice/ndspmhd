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
 
 use kernels,      only:interpolate_kernels,radkern2
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

 integer, parameter :: maxsweeps = 1000
 integer :: i,j,n,nsweeps
 integer :: icell,iprev,nneigh,nneighi
 integer, dimension(npart) :: listneigh ! neighbour list
 real :: rij,rij2
 real :: hi1,hi21,etaij,etai,etaj
 real :: hfacwabi,rho1j,sumdenom
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,Bi,dB,sumB,sumBin,lhs,rhs
 real :: q2i,grkerni,wabi,dB2,dBmax,term,fac,errmax
 real :: grkernalti,grgrkernalti,grgrkerni,dti
 real, dimension(ndimV,ntotal) :: Bfieldnew
 logical :: converged
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
 nsweeps = 0
 converged = .false.
!
!--fac is the explicit/implicit factor
!  fac = 0 gives fully explicit (forwards Euler)
!  fac = 1 gives fully implicit (backwards Euler)
!  fac = 1/2 gives mix (Crank-Nicolson)
 
 fac = 0.5
 if (abs(dt) > 0.) then
    dti = dt
 else
    !--if dt=0, return explicit derivative
    dti = 1.
    fac = 0.
 endif
 
 sweeps: do while (.not.converged .and. nsweeps.lt.maxsweeps)
 
 nsweeps = nsweeps + 1
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
 errmax = 0.
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
       etai = etafunc(x(1,i),etamhd)
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
                etaj = etafunc(x(1,j),etamhd)
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                !  (using hi)
                !call interpolate_kernel(q2i,wabi,grkerni)
                call interpolate_kernels(q2i,wabi,grkerni,grkernalti,grgrkernalti)
!
!--calculate nabla^2 B
!
                grkerni = grkerni*hfacwabi*hi1*gradh(i)
                grkernalti = grkernalti*hfacwabi*hi1

                !--usual 2nd deriv
                grgrkerni = -2.*grkerni/rij
                
                !--uncomment following lines to use alt kernel
                !grgrkerni = -2.*grkernalti/rij
                
                !--uncomment following lines to use 2nd deriv of alternative kernel
                !grgrkerni = grgrkernalti*hfacwabi*hi1*hi1
                
                !--constant resistivity
                !etaij = etamhd
                !--arithmetic average
                !etaij = 0.5*(etai + etaj)
                !--geometric average
                etaij = 2.*etai*etaj/(etai + etaj)

                term = etaij*pmass(j)*rho1j*grgrkerni
                !term = etamhd*pmass(j)*rho1j*grgrkerni
                sumdenom = sumdenom + term
                sumB(:) = sumB(:) + term*Bfieldnew(:,j)
                sumBin(:) = sumBin(:) + term*Bfield(:,j)
                nneighi = nneighi + 1

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours

       Bi(:) = Bfieldnew(:,i)
       Bfieldnew(:,i) = (Bfield(:,i)*(1. - (1.-fac)*dti*sumdenom) &
                         + fac*dti*sumB(:) + (1.-fac)*dti*sumBin(:))/&
                        (1. + fac*dti*sumdenom)
       !print*,' particle ',i,' numneigh = ',numneigh(i),nneighi
       
       dB = Bfieldnew(:,i) - Bi(:)
       dB2 = dot_product(dB,dB)
       dBmax = max(dBmax,dB2)
       
       !
       !--check LHS=RHS convergence
       !
       if (dt.gt.0.) then
          lhs = Bfieldnew(:,i) ! - Bfield(:,i))/dt
          rhs = Bfield(:,i) + dt*((1.-fac)*(sumBin(:) - sumdenom*Bfield(:,i)) &
                   + fac*(sumB(:) - sumdenom*Bfieldnew(:,i)))
          errmax = max(maxval(abs(rhs-lhs)),errmax)
       endif

       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells
 
 converged = (dBmax < 1.e-12 .and. errmax < 1.e-12)
 !print*,nsweeps,' errmax = ',errmax
 !print*,' end of sweep ',nsweeps,' dBmax = ',dBmax,' dt =  ',dt
 if (fac < 1. .and. nsweeps.eq.maxsweeps .and. .not.converged) then
    print*,' NOT CONVERGED setting fac = 1 and retrying: press enter to continue'
    stop
    !read*
    do i=1,ntotal
       Bfieldnew(:,i) = Bfield(:,i)
    enddo
    fac = 1.
    nsweeps = 0
 endif

 enddo sweeps
 
 if (converged) then
    print*,' converged in ',nsweeps,' sweeps'
 else
    print*,' ERROR: resistivity not converged, max error = ',dBmax
    !stop
 endif
 !read*

!
!--what we actually return from this routine is an extra contribution to
!  dBevoldt, as if it had been computed with the B^(n+1) in it.
!
    !dBevoldt(:,:) = 0.
 if (dti.gt.0.) then
    do i=1,npart
       dBevoldt(:,i) = dBevoldt(:,i) + (Bfieldnew(:,i) - Bfield(:,i))/dti
    enddo
 endif

 return
end subroutine Bdiffusion

!
!--diffusion parameter as a function of x
!
real function etafunc(x,eta)
 implicit none
 real, parameter :: pi = 3.1415926536
 real, intent(in) :: x,eta
 
! if (x > 0.25 .and. x < 0.75) then
!    etafunc = 10000.*eta
! else
!    etafunc = 0.01*eta
! endif
 etafunc = eta*cos(2.*pi*x)**2
 
end function etafunc

end module resistivity
