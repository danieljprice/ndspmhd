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
! Computes terms required for Quantum SPH
!
! Input:
!    iderivtype - type of second derivative operator (see below)
!    npart - number of particles
!    x     - positions
!    pmass - particle masses
!    rho   - density
!    hh    - smoothing length
!
! Output:
!    P_Q(:,:,npart) - quantum pressure tensor
!
! iderivtype = 1 (default): Brookshaw/Espanol-Revenga style derivatives
!
!------------------------------------------------------------------------
module get_quantum
 implicit none
 
contains

subroutine get_quantum_pressure(iderivtype,npart,x,pmass,rho,hh,P_Q)
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
 real, dimension(ndim,ndim,idim), intent(out) :: P_Q

 integer :: i,j,k,l,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: rho1i,rho1j,rij1,rhoi,rhoj,drho
 real, dimension(ntotal) :: h1
 real, dimension(ndim,ntotal) :: gradrho
 real, dimension(ndim,ndim,ntotal) :: gradgradrho
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_quantum' 
!
!--initialise quantities
!
 listneigh = 0
 dr(:) = 0.
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
 gradrho(:,:) = 0.
 gradgradrho(:,:,:) = 0.
 P_Q(:,:,:) = 0.
 select case(iderivtype)
 case(2)
    print "(a)",' computing Quantum pressure with different method'
 case default
    print "(a)",' computing Quantum pressure with default method' 
 end select
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
                call interpolate_kernel_curl(q2i,grkerni,grgrkerni)
                !  (using hj)
                call interpolate_kernel_curl(q2j,grkernj,grgrkernj)
!
!--calculate 2nd derivatives of rho
!
                rhoj  = rho(j)
                rho1j = 1./rhoj
                drho = rhoi - rhoj

                select case(iderivtype)
                case(2)  ! "alternative method"
                   ! insert favourite method here
                case default  ! "default method"
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1
!
!--1st derivatives of density
!
                   gradrho(:,i) = gradrho(:,i) + pmass(j)*dr(1:ndim)*grkerni
                   gradrho(:,j) = gradrho(:,j) - pmassi*dr(1:ndim)*grkernj
!
!--2nd derivatives of density
!
                   do k=1,ndim
                      do l=1,ndim
                         gradgradrho(l,k,i) = gradgradrho(l,k,i) + &
                            pmass(j)/rho(j) * (rho(i) - rho(j))*((ndim+2)*dr(l)*dr(k) - delta_fn(l,k,ndim))*grkerni*rij1
                         gradgradrho(l,k,j) = gradgradrho(l,k,j) + &
                           pmassi/rho(i) * (rho(j) - rho(i))*((ndim+2)*dr(l)*dr(k) - delta_fn(l,k,ndim))*grkernj*rij1
                      enddo
                   enddo
                end select
             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells
 
 do n=1,npart
    do j=1,ndim
       do i=1,ndim
       	  !P_Q(i,j,n) = gradrho(i, n)
          !P_Q(i,j,n) = gradgradrho(i,j,n)
          P_Q(i,j,n) = - 0.25*(gradgradrho(i,j,n) - gradrho(i,n)*gradrho(j,n)/rho(n))
       enddo
    enddo
 enddo
 
 !print*,' P_Q = ',P_Q(:,:,1:10)

 return
end subroutine get_quantum_pressure

integer function delta_fn(i,j,ndim)
 integer, intent(in) :: i,j,ndim
 
 if (i==j) then
    delta_fn = ndim
 else
    delta_fn = 0
 endif
 
end function delta_fn
     
end module get_quantum
