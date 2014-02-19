!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
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

module getdivB
 use get_neighbour_lists
 implicit none
 
contains

!!------------------------------------------------------------------------
!! Computes an SPH estimate of div B
!! This version computes div B on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

subroutine get_divB(divBonrho,ntot)
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 
 use hterms, only:gradh
 use kernels, only:interpolate_kernel,radkern2
 use linklist
 use options, only:ikernav
 use part
 use setup_params
!
!--define local variables
!
 implicit none
 integer, intent(in) :: ntot
 real, dimension(:), intent(out) :: divBonrho
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
 real :: rho21i, rho21j, projdb
 real, dimension(ndim) :: dx
 real, dimension(ndimv) :: dr
!
!  (kernel quantities)
!
 real :: q2,q2i,q2j      
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine get_divB, ntot=',ntot  
!
!--initialise quantities
!
 !!!print*,'ntot = ',ntot, 'size= ',size(divBonrho)
 do i=1,ntot
    divBonrho(i) = 0.
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
    call get_neighbour_list(icell,listneigh,nneigh)
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
       rho21i = 1./rho(i)**2
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: DO n = idone+1,nneigh
          j = listneigh(n)
          if ((j.NE.i).AND..NOT.(j.GT.npart .AND. i.GT.npart)) then
             ! don't count particle with itself       
!!            
             if (j.lt.0 .or. j.gt.ntotal) then
                print*,i,j
                print*,listneigh(1:nneigh)
                print*,'n=',n,' listneigh(n) = ',listneigh(n)
             endif
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
             
             rij2 = DOT_PRODUCT(dx,dx)
             rij = SQRT(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2     
             dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             if (ndimV.gt.ndim) dr(ndim+1:ndimV) = 0. 

             !          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2),dx,hi,hj
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if (((q2i.lt.radkern2).OR.(q2j.lt.radkern2))  &
                  .and. .not.(i.gt.npart.and.j.gt.npart)) then
!     
!--interpolate from kernel table          
!  (use either average h or average kernel gradient)
!
                !       PRINT*,' neighbour,r/h,dx,hi,hj ',i,j,SQRT(q2),dx,hi,hj
                if (ikernav.EQ.1) then          
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

                if (ikernav.NE.3) then
                   grkerni = grkern
                   grkernj = grkern
                endif         
!
!--calculate div B
!
                !term = DOT_PRODUCT((Bfield(:,i)*rho21i*grkerni*gradh(i) &
                !                  + Bfield(:,j)*rho21j*grkernj*gradh(j)),dr)
                projdB = DOT_PRODUCT(Bfield(:,i)-Bfield(:,j),dr)

                divBonrho(i) = divBonrho(i) - pmass(j)*projdB*grkerni
                divBonrho(j) = divBonrho(j) - pmass(i)*projdB*grkernj            

                !divBonrho(i) = divBonrho(i) + pmass(j)*term
                !divBonrho(j) = divBonrho(j) - pmass(i)*term
                !      else
                !         PRINT*,' r/h > 2 '      
                
             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.NE.-1) i = ll(i)          ! possibly should be only if (iprev.NE.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 if (ikernav.EQ.3) then
    do i=1,ntotal
       divBonrho(i) = gradh(i)*divBonrho(i)/rho(i)**2
    enddo
 else
    do i=1,ntotal
       divBonrho(i) = divBonrho(i)/rho(i)**2
    enddo
 endif

 return
end subroutine get_divb

end module getdivB   
