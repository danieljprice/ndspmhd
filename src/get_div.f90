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

!!------------------------------------------------------------------------
!! Computes an SPH estimate of div B given a vector B
!!
!! This version computes div B on all particles
!! and therefore only does each pairwise interaction once
!!
!! particle quantities are explicitly passed rather than through the 
!! global modules
!!
!! assumes there has been a previous call to density
!!------------------------------------------------------------------------

subroutine get_divBgradpsi(divb,bin,x,hh,pmass,rho,npart,ntot)
 use dimen_mhd, only:ndim, ndimV
 use debug, only:trace
 use loguns, only:iprint
 
 use bound, only:ireal
 use kernels, only:interpolate_kernel
 use linklist
 use options, only:ikernav
 use get_neighbour_lists
 implicit none
 integer, intent(in) :: npart,ntot
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: hh,pmass,rho
 real, dimension(ndimv,ntot), intent(in) :: bin
 real, dimension(ntot), intent(out) :: divb
!
!--define local variables
!
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
 real :: pmassi,pmassj,projdb
 real, dimension(ndim) :: dx
 real, dimension(ndimv) :: dr
!
!  (kernel quantities)
!
 real :: q2,q2i,q2j      
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj,weight
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_divbgradpsi, npart,ntot=',npart,ntot 
!
!--initialise quantities
!
 listneigh = 0
 divb = 0.
 rho = 0.
!
!--loop over all the link-list cells
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
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1) ! loop over home cell particles

       !print*,'doing particle ',i,nneigh,' neighbours',hh(i),bin(:,i)
       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hi2 = hi*hi
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if (.not.(j.gt.npart .and. i.gt.npart)) then
             ! do count particle with itself       
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
             pmassj = pmass(j)
             
             rij2 = dot_product(dx,dx)
             rij = sqrt(rij2)
             q2 = rij2/h2
             q2i = rij2/hi2
             q2j = rij2/hj2
             dr = 0.
             if (j.ne.i) then
                dr(1:ndim) = dx(1:ndim)/rij  ! unit vector
             endif
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
                !print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                if (ikernav.eq.1) then          
                   call interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hj1
                   grkerni = grkern
                   grkernj = grkern
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
                   if (ikernav.eq.2) then
                      wabi = wab
                      wabj = wab
                      grkerni = grkern
                      grkernj = grkern
                   endif           
                endif
!
!--calculate div b and grad psi
!
                projdb = dot_product(bin(:,i)-bin(:,j),dr)

                divb(i) = divb(i) - pmassj*projdb*grkerni
                divb(j) = divb(j) - pmassi*projdb*grkernj            
                                
             endif

         endif! .not. j>npart .and. i>npart   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

!
!--do divisions by rho
!  note that grad psi returns the correct term appropriate to evolving either
!  b or b/rho (ie. grad psi or grad psi / rho respectively). this should be
!  the same in rates.
!
 do i=1,npart
    if (rho(i).ge.0.) then
       divb(i) = divb(i)/(gradh(i)*rho(i))
    else
       write(*,*) 'error: get_divbgradpsi: rho < 0.'
    endif
    !print*,'divb, rho = ',divb(i),rho(i),bin(:,i)
    !if (mod(i,20).eq.0) read*
 enddo
 
 do i=npart+1,ntot
    j = ireal(i)
    divb(i) = divb(j)
 enddo
! read*

 return
end subroutine get_divbgradpsi
      
