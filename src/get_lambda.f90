module getlambda
 use get_neighbour_lists
 implicit none
 
contains

!!------------------------------------------------------------------------
!! Get the value of lambda from a linear combination of two kernels
!!------------------------------------------------------------------------

subroutine get_lambda(lambda,ntot)
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 
 use hterms, only:gradh
 use kernels, only:interpolate_kernels2,radkern2
 use linklist
 use options, only:ikernav
 use part
 use setup_params
!
!--define local variables
!
 implicit none
 integer, intent(in) :: ntot
 real, dimension(:), intent(out) :: lambda
 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(ntot) :: listneigh ! neighbour list
 integer :: idone
!
!  (particle properties - local copies)
!      
 real :: rij,rij2
 real :: hi,hi1,hj,hj1,h2,hi2,hj2
 real :: hfacwabi,hfacwabj
 real :: rho1i,rho1j,rho21i, rho21j
 real :: rhoalt1i,rhoalt1j,rhoalt21i,rhoalt21j
 real, dimension(ndim) :: dx
 real, dimension(ndimv) :: dr
!
!  (kernel quantities)
!
 real :: q2,q2i,q2j      
 real :: wabi,wabj
 real :: wabalti,wabaltj
 real :: grkerni,grkernj,grkernalti,grkernaltj
 real :: term,termalt
 real, dimension(ndimV) :: gradunityi,gradunityalti
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine get_lambda, ntot=',ntot  
!
!--initialise quantities
!
 !!!print*,'ntot = ',ntot, 'size= ',size(divBonrho)
 do i=1,ntot
    lambda(i) = 0.
    rhoalt(i) = pmass(i)*rhoalt(i)
 enddo
 listneigh = 0
!
!--Loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop          ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list_partial(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1     ! note density summation includes current particle
    if (i.NE.-1) iprev = i

    loop_over_cell_particles: do while (i.NE.-1)          ! loop over home cell particles

       !       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       rho1i  = 1./rho(i)
       rho21i = rho1i*rho1i
       rhoalt1i  = 1./rhoalt(i)
       rhoalt21i = rhoalt1i*rhoalt1i 
       gradunityi(:) = 0.
       gradunityalti(:) = 0.
       if (i.gt.npart) then
          i = ll(i)
          cycle loop_over_cell_particles
       endif
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = 1,nneigh
          j = listneigh(n)
          if (j.lt.0 .or. j.gt.ntotal) then
             print*,i,j
             print*,listneigh(1:nneigh)
             print*,'n=',n,' listneigh(n) = ',listneigh(n)
          endif
          dx(:) = x(:,i) - x(:,j)
          hj = hh(j)
          hj1 = 1./hj
          hj2 = hj*hj
          hfacwabj = hj1**ndim
          rho1j  = 1./rho(j)
          rho21j = rho1j*rho1j
          rhoalt1j  = 1./rhoalt(j)
          rhoalt21j = rhoalt1j*rhoalt1j

          rij2 = DOT_PRODUCT(dx,dx)
          q2 = rij2/h2
          q2i = rij2/hi2
          q2j = rij2/hj2
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
          if (((q2i.lt.radkern2).OR.(q2j.lt.radkern2))  &
               .and. .not.(i.gt.npart.and.j.gt.npart)) then

             if (rij2.gt.0.) then
                rij = SQRT(rij2)
                dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             else
                dr(:) = 0.
             endif
             if (ndimV.gt.ndim) dr(ndim+1:ndimV) = 0.
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
             !  (using hi)
             call interpolate_kernels2(q2i,wabi,wabalti,grkerni,grkernalti)
             wabi    = wabi*hfacwabi
             wabalti    = wabalti*hfacwabi
             grkerni    = grkerni*hfacwabi*hi1
             grkernalti = grkernalti*hfacwabi*hi1
             !  (using hj)
             call interpolate_kernels2(q2j,wabj,wabaltj,grkernj,grkernaltj)
             wabj = wabj*hfacwabj
             wabaltj = wabaltj*hfacwabj
             grkernj = grkernj*hfacwabj*hj1
             grkernaltj = grkernaltj*hfacwabj*hj1
!
!--calculate term
!
             !term    = pmass(j)*wabi*rho1j
             !termalt = pmass(j)*wabalti*rhoalt1j
             term    = pmass(j)*(rho21i*grkerni + rho21j*grkernj)
             termalt = pmass(j)*(rhoalt21i*grkernalti + rhoalt21j*grkernaltj)

             gradunityi(:)    = gradunityi(:)    + term*dr(:)
             gradunityalti(:) = gradunityalti(:) + termalt*dr(:) 
          endif
       enddo loop_over_neighbours
!
!--compute lambda for this particle
!
       lambda(i) = gradunityi(1)
       if (i.eq.50) then
          print*,i,'h = ',hh(i)
          print*,i,'lambda = ',rho(i)*gradunityi(1),rhoalt(i)*gradunityalti(1),rho(i),rhoalt(i),gradh(i)
       endif
       iprev = i
       if (iprev.NE.-1) i = ll(i)          ! possibly should be only if (iprev.NE.-1)
    enddo loop_over_cell_particles

 enddo loop_over_cells

 return
end subroutine get_lambda

end module getlambda
