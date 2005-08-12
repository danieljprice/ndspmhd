!!------------------------------------------------------------------------
!! Computes an SPH estimate of the curl of a vector quantity
!! This version computes the curl for all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

subroutine get_curl(curlbonrho,ntot)
 use dimen_mhd, only:ndim,ndimv
 use debug, only:trace
 use loguns, only:iprint
 
 use bound
 use kernels, only:interpolate_kernel,radkern2
 use linklist
 use options, only:ikernav,imhd
 use part
 use setup_params
 use get_neighbour_lists
!
!--define local variables
!
 implicit none
 integer, intent(in) :: ntot
 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(ntot) :: listneigh ! neighbour list
 integer :: idone
!
!  (particle properties - local copies)
!      
 real :: rij,rij2
 real :: hi,hi1,hav,hav1,hj,hj1,h2,hi2,hj2
 real :: hfacwab,hfacwabi,hfacwabj
 real :: rho21i, rho21j, pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimv) :: dr, db, bi, bj, curlbi
 real, dimension(ndimv,ntot), intent(out) :: curlbonrho
!
!  (kernel quantities)
!
 real :: q2,q2i,q2j      
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_curl, ntot=',ntot  
!
!--initialise quantities
!
 listneigh = 0
 curlbonrho = 0.
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
       rho21i = 1./rho(i)**2
       pmassi = pmass(i)
       if (imhd.ge.11) then      ! if mag field variable is b
          bi(:) = bevol(:,i)
       elseif (imhd.gt.0) then      ! if mag field variable is b/rho
          bi(:) = bfield(:,i)
       endif
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
!
!--calculate averages of smoothing length if using this averaging
!                
             hav = 0.5*(hi + hj)
             hav1 = 1./hav
             h2 = hav*hav
             hfacwab = hav1**ndim
             hfacwabj = hj1**ndim
             rho21j = 1./rho(j)**2
             
             rij2 = dot_product(dx,dx)
             rij = sqrt(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2     
             dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             if (ndimv.gt.ndim) dr(ndim+1:ndimv) = 0. 

             !          print*,' neighbour,r/h,dx,hi,hj ',j,sqrt(q2),dx,hi,hj
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if (((q2i.lt.radkern2).or.(q2j.lt.radkern2))  &
                  .and. .not.(i.gt.npart.and.j.gt.npart)) then
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                if (ikernav.eq.1) then          
                   call interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hj1
                else
                   !  (using hi)
                   call interpolate_kernel(q2i,wabi,grkerni)
                   wabi = wabi*hfacwabi
                   grkerni = grkerni*hfacwabi*hi1
                   !  (using hj)
                   call interpolate_kernel(q2j,wabj,grkernj)
                   wabj = wabj*hfacwabj
                   grkernj = grkernj*hfacwabj*hj1
                   !  (calculate average)            
                   wab = 0.5*(wabi + wabj)                  
                   grkern = 0.5*(grkerni + grkernj)            
                endif

                if (imhd.ge.11) then      ! if b is mag field variable
                   bj(:) = bevol(:,j)
                elseif (imhd.ne.0) then      ! if b/rho is mag field variable
                   bj(:) = bfield(:,j)                                
                endif
!
!--calculate div b
!
                db = bi(:) - bj(:)
                if (ndimv.eq.3) then
                   curlbi(1) = db(2)*dr(3) - db(3)*dr(2)
                   curlbi(2) = db(3)*dr(1) - db(1)*dr(3)
                   curlbi(3) = db(1)*dr(2) - db(2)*dr(1)
                elseif (ndimv.eq.2) then  ! just jz in 2d
                   curlbi(1) = db(1)*dr(2) - db(2)*dr(1)
                   curlbi(2) = 0.
                endif
                !
                !--compute rho * current density j
                !
                curlbonrho(:,i) = curlbonrho(:,i) - pmass(j)*curlbi(:)*grkern
                curlbonrho(:,j) = curlbonrho(:,j) - pmassi*curlbi(:)*grkern
                !      else
                !         print*,' r/h > 2 '      
                
             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,ntot
    curlbonrho(:,i) = curlbonrho(:,i)/rho(i)**2
 enddo

 return
end subroutine get_curl
      
