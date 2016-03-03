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
! Computes an SPH estimate of the divergence of a vector quantity
! This version computes the divergence for all particles
! and therefore only does each pairwise interaction once
!
! Input:
!    idivtype - type of divergence operator (see below)
!    npart - number of particles
!    x     - positions
!    pmass - particle masses
!    rho   - density
!    hh    - smoothing length
!    a     - scalar quantity to take divergence of
!
! Output:
!    diva - divergence of a
!
! idivtype = 1 (default): mass-weighted differenced divergence operator:
! 
!    -1/(rho_i*Omega_i) \sum m_j (A_i - A_j) .\nabla W_ij (h_i)
!
! idivtype = 2: mass-weighted symmetric divergence operator:
! 
!    rho_i \sum m_j ((A_i .\nabla W_ij(h_i)) / (rho_i^2 Omega_i) +
!                   ((A_j .\nabla W_ij(h_j)) / (rho_j^2 Omega_j)
!    
! idivtype = 3: differenced divergence operator with constant weights:
! 
!    -m_i/(rho_i*h_i**3) \sum (A_i - A_j) \nabla W_ij (h_i)
!
! idivtype = 4: differenced divergence operator:
! 
!    -rho_i \sum m_j/rho_j**2 (A_i - A_j) \nabla W_ij (h_i)
!
!------------------------------------------------------------------------
module getdiv
 implicit none
 
contains

subroutine get_divergence(idivtype,npart,x,pmass,rho,hh,a,diva)
 use dimen_mhd,           only:ndim,ndimV
 use debug,               only:trace
 use loguns,              only:iprint
 use kernels,             only:interpolate_kernel,radkern2
 use linklist,            only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms,              only:gradh
 use setup_params,        only:hfact
 use part,                only:ntotal
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart,idivtype
 real, dimension(:,:), intent(in)  :: x
 real, dimension(:),   intent(in)  :: pmass,rho,hh
 real, dimension(:,:), intent(in)  :: a
 real, dimension(:),   intent(out) :: diva
 real :: weight

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi,divterm,divai
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,ai
 real :: q2i,q2j,wi,wj,grkerni,grkernj
 real :: rho21i,rho21gradhi,termi,termj
 real, dimension(ntotal) :: h1
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_div' 
!
!--initialise quantities
!
 listneigh = 0
 diva = 0.
 dr(:) = 0.
 weight = 1./hfact**ndim
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
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
       !!hi = hh(i)
       hi1 = h1(i) !!1./hi
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       pmassi = pmass(i)
       ai     = a(:,i)
       rho21i = 1./rho(i)**2
       rho21gradhi = rho21i*gradh(i)
       divai  = 0.
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             !!hj = hh(j)
             hj1 = h1(j) !!1./hj
             
             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
             q2j = rij2*hj1*hj1
!     
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then

                hfacwabj = hj1**ndim
                rij = sqrt(rij2)
                dr(1:ndim) = dx(1:ndim)/(rij + tiny(rij))  ! unit vector
!     
!--interpolate from kernel table          
!  (use either average h or average kernel divergence)
!
                !  (using hi)
                call interpolate_kernel(q2i,wi,grkerni)
                !  (using hj)
                call interpolate_kernel(q2j,wj,grkernj)
!
!--calculate divergence of A
!
                select case(idivtype)
                case(2)  ! symmetric divergence operator
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1*gradh(j)

                   termi = dot_product(ai(1:ndim),dr(1:ndim))
                   termj = dot_product(a(:,j),dr(1:ndim))
                   divterm   = termi*rho21gradhi*grkerni + termj/rho(j)**2*grkernj
                   divai   = divai   + pmass(j)*divterm
                   diva(j) = diva(j) - pmassi*divterm
                
                case(3)  ! (dB using weights) -- multiply by weight below
                   grkerni = grkerni*hi1
                   grkernj = grkernj*hj1

                   divterm = dot_product(ai(1:ndim) - a(1:ndim,j),dr(1:ndim))
                   divai   = divai   + divterm*grkerni
                   diva(j) = diva(j) + divterm*grkernj

                case(4)  ! (dB, m_j/rho_j**2) -- multiply by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   divterm = dot_product(ai(1:ndim) - a(1:ndim,j),dr(1:ndim))
                   
                   divai   = divai   + pmass(j)/rho(j)**2*divterm*grkerni
                   diva(j) = diva(j) + pmassi*rho21i*divterm*grkernj

                case default  ! default differenced divergence -- divide by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   divterm = dot_product(ai(1:ndim) - a(1:ndim,j),dr(1:ndim))

                   divai   = divai   + pmass(j)*divterm*grkerni
                   diva(j) = diva(j) + pmassi*divterm*grkernj
                end select

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       diva(i) = diva(i) + divai
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    select case(idivtype)
    case(4)
       diva(i) = -rho(i)*diva(i)
    case(3)
       diva(i) = -weight*diva(i) !!*gradh(i)
    case(2)
       diva(i) = rho(i)*diva(i)
    case default
       diva(i) = -diva(i)*gradh(i)/rho(i)
    end select
 enddo

 return
end subroutine get_divergence
      
end module getdiv
