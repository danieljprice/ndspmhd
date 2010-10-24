!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2010 Daniel Price                                                   !
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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Smoothes out initial conditions using the SPH summation, that is the  !!
!!   quantity x_a on particle a is recalculated by summing over the       !!
!!   particle neighbours according to:                                    !!
!!                                                                        !!
!!   x_smooth(r_a) = sum_b m_b (x_b / rho_b) W(r-r_b,h_ab)                 !!
!!                                                                        !!
!!   where W is the SPH smoothing kernel                                   !!
!!                                                                          !!
!!  Needs to be preceded by a call to linklist to setup neighbour lists   !!
!!  The density should also be previously known                                  !!
!!------------------------------------------------------------------------!!
module smooth
 implicit none
 public :: smooth_variable
 
 private
 
contains

subroutine smooth_variable(xinput,xsmooth,x,pmass,hh,rho)
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use kernels, only:interpolate_kernel,radkern2
 use linklist, only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists
 use part, only:npart,ntotal
 implicit none
 real, dimension(:), intent(in) :: xinput
 real, dimension(:,:), intent(in) :: x
 real, dimension(:), intent(in) :: pmass,hh,rho
 real, dimension(:), intent(out) :: xsmooth
!
!--define local variables
!      
 integer :: i,j,n
 integer, dimension(ntotal) :: listneigh ! neighbour list
 integer :: idone,iprev,nneigh,icell
 real, dimension(ndim) :: dx
 real :: termi,termj,hi,hi21
 real :: q2i,q2j,rij2,weight,wabi,wabj,grkerni,grkernj
!
!--allow for tracing flow
!
 write(iprint,*) 'applying smoothing...'
!
!--initialise quantities
!
 listneigh = 0
 xsmooth(:) = 0.
!
!--check size of input array is consistent with rho
!
 if (size(xinput).ne.size(rho)) then 
    write(iprint,*)  &
      'error: smooth: input array and rho have different dimensions'
    call quit   
 endif   
!
!--loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop                ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
    call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
    i = ifirstincell(icell)
    idone = -1        ! note density summation includes current particle
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1)                ! loop over home cell particles

!       print*,'doing particle ',i,nneigh,' neighbours',pmass(i)
       idone = idone + 1
       hi = hh(i)
       termi = pmass(i)*xinput(i)/(rho(i)*hi**ndim)
       hi21 = 1./hi**2
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if (.not.(i.gt.npart.and.j.gt.npart)) then
             dx(:) = x(:,i) - x(:,j)
             termj = pmass(j)*xinput(j)/(rho(j)*hh(j)**ndim)
             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
             q2j = rij2/hh(j)**2
!        
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             weight = 1.0
             if (j.eq.i) weight = 0.5
             if (q2i.lt.radkern2) then
!        
!--interpolate from kernel table (using hi and hj)
!
                call interpolate_kernel(q2i,wabi,grkerni)
                xsmooth(i) = xsmooth(i) + termj*wabi*weight
             endif
             if (q2j.lt.radkern2) then
!        
!--interpolate from kernel table (using hi and hj)
!
                call interpolate_kernel(q2j,wabj,grkernj)
                   xsmooth(j) = xsmooth(j) + termi*wabj*weight
             endif
         endif
            
       enddo loop_over_neighbours
            
       iprev = i
       if (iprev.ne.-1) i = ll(i)   ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
            
 enddo loop_over_cells
!
!--update ghosts
! 
 if (ntotal.gt.npart) then
    do i=npart+1,ntotal
       j = ireal(i)
       xsmooth(i) = xsmooth(j)
    enddo
 endif   

 return
end subroutine smooth_variable

end module smooth
