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

!------------------------------------------------------------------------
! Computes an SPH estimate of the curl of a vector quantity
! This version computes the curl for all particles
! and therefore only does each pairwise interaction once
!
! Input:
!    icurltype - type of curl operator (see below)
!    npart - number of particles
!    x     - positions
!    pmass - particle masses
!    rho   - density
!    hh    - smoothing length
!    Bvec  - vector quantity to take curl of
!
! Output:
!    curlB - curl of Bvec
!    (optional) curlBgradh - term required for vector potential stuff
!
! icurltype = 1 (default): mass-weighted differenced curl operator:
! 
!    1/(rho_i*Omega_i) \sum m_j (B_i - B_j) x \nabla W_ij (h_i)
!
! icurltype = 2: mass-weighted symmetric curl operator:
! 
!    -rho_i \sum m_j ((B_i x \nabla W_ij(h_i)) / (rho_i^2 Omega_i) +
!                    ((B_j x \nabla W_ij(h_j)) / (rho_j^2 Omega_j)
!    
! icurltype = 3: differenced curl operator with constant weights:
! 
!    m_i/(rho_i*h_i**3) \sum (B_i - B_j) x \nabla W_ij (h_i)
!
! icurltype = 4: differenced curl operator:
! 
!    rho_i \sum m_j/rho_j**2 (B_i - B_j) x \nabla W_ij (h_i)
!
!------------------------------------------------------------------------
module getcurl
 implicit none
 
contains

subroutine get_curl(icurltype,npart,x,pmass,rho,hh,Bvec,curlB,curlBgradh)
 use dimen_mhd,           only:ndim,ndimV,idim
 use debug,               only:trace
 use loguns,              only:iprint
 use kernels,             only:interpolate_kernel_curl,radkern2
 use linklist,            only:ll,ifirstincell,ncellsloop
 use get_neighbour_lists, only:get_neighbour_list
 use hterms,              only:gradh
 use setup_params,        only:hfact
 use part,                only:ntotal
 use utils,               only:cross_product3D
!
!--define local variables
!
 implicit none
 integer, intent(in) :: npart,icurltype
 real, dimension(ndim,idim), intent(in) :: x
 real, dimension(idim), intent(in) :: pmass,rho,hh
 real, dimension(ndimV,idim), intent(in) :: Bvec
 real, dimension(ndimV,idim), intent(out) :: curlB
 real, dimension(ndimV,idim), intent(out), optional :: curlBgradh
 real :: weight

 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, dimension(npart) :: listneigh ! neighbour list
 integer :: idone
 real :: rij,rij2
 real :: hi1,hj1,hi21
 real :: hfacwabi,hfacwabj
 real :: pmassi
 real, dimension(ndim) :: dx
 real, dimension(ndimV) :: dr,dB,Bi,curlBi,curlBterm,curlBtermi,curlBtermj,curlBgradhi
 real :: q2i,q2j,grkerni,grkernj,grgrkerni,grgrkernj
 real :: dgradwdhi,dgradwdhj,rho21i,rho21gradhi
 real, dimension(ntotal) :: h1
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine get_curl' 
!
!--initialise quantities
!
 listneigh = 0
 curlB = 0.
 dr(:) = 0.
 weight = 1./hfact**ndim
 do i=1,ntotal
    h1(i) = 1./hh(i)
 enddo
 if (present(curlBgradh)) curlBgradh(:,:) = 0.
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
       Bi(:) = Bvec(:,i)
       rho21i = 1./rho(i)**2
       rho21gradhi = rho21i*gradh(i)
       curlBi(:) = 0.
       curlBgradhi(:) = 0.
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
             !          print*,' neighbour,r/h,dx,hi,hj ',j,sqrt(q2),dx,hi,hj
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
!  (use either average h or average kernel gradient)
!
                !       print*,' neighbour,r/h,dx,hi,hj ',i,j,sqrt(q2),dx,hi,hj
                !  (using hi)
                call interpolate_kernel_curl(q2i,grkerni,grgrkerni)
                !  (using hj)
                call interpolate_kernel_curl(q2j,grkernj,grgrkernj)
!
!--calculate curl of Bvec (NB dB is 3-dimensional, dr is also but zero in parts)
!
                select case(icurltype)
                case(2)  ! symmetric curl for vector potential current
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1*gradh(j)

                   call cross_product3D(Bi(:),dr,curlBtermi)
                   call cross_product3D(Bvec(:,j),dr,curlBtermj)
                   curlBterm(:) = (curlBtermi(:)*rho21gradhi*grkerni + curlBtermj(:)/rho(j)**2*grkernj)
                   curlBi(:) = curlBi(:) + pmass(j)*curlBterm(:)
                   curlB(:,j) = curlB(:,j) - pmassi*curlBterm(:)
                
                case(3)  ! (dB using weights) -- multiply by weight below
                   grkerni = grkerni*hi1
                   grkernj = grkernj*hj1

                   dB(1:ndimV) = Bi(1:ndimV) - Bvec(1:ndimV,j)
                   call cross_product3D(dB,dr,curlBterm)
                   curlBi(:) = curlBi(:) + curlBterm(:)*grkerni
                   curlB(:,j) = curlB(:,j) + curlBterm(:)*grkernj

                case(4)  ! (dB, m_j/rho_j**2) -- multiply by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   dB(1:ndimV) = Bi(1:ndimV) - Bvec(1:ndimV,j)
                   call cross_product3D(dB,dr,curlBterm)
                   
                   curlBi(:) = curlBi(:) + pmass(j)/rho(j)**2*curlBterm(:)*grkerni
                   curlB(:,j) = curlB(:,j) + pmassi*rho21i*curlBterm(:)*grkernj

                case default  ! default curl (dB, m_j/rho_i with gradh) -- divide by rho(i) below
                   grkerni = grkerni*hfacwabi*hi1
                   grkernj = grkernj*hfacwabj*hj1

                   dB(1:ndimV) = Bi(1:ndimV) - Bvec(1:ndimV,j)
                   call cross_product3D(dB,dr,curlBterm)

                   curlBi(:) = curlBi(:) + pmass(j)*curlBterm(:)*grkerni
                   curlB(:,j) = curlB(:,j) + pmassi*curlBterm(:)*grkernj

                   if (present(curlBgradh)) then
                      dgradwdhi = -(ndim+1.)*hi1*grkerni - rij*hi1**3*grgrkerni*hfacwabi
                      dgradwdhj = -(ndim+1.)*hj1*grkernj - rij*hj1**3*grgrkernj*hfacwabj

                      curlBgradhi(:) = curlBgradhi(:) + pmass(j)*curlBterm(:)*dgradwdhi
                      curlBgradh(:,j) = curlBgradh(:,j) + pmassi*curlBterm(:)*dgradwdhj
                   endif
                end select

             endif
          endif! j .ne. i   
       enddo loop_over_neighbours
       
       curlB(:,i) = curlB(:,i) + curlBi(:)
       if (present(curlBgradh)) then
          curlBgradh(:,i) = curlBgradh(:,i) + curlBgradhi(:)
       endif
       iprev = i
       if (iprev.ne.-1) i = ll(i)          ! possibly should be only if (iprev.ne.-1)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=1,npart
    select case(icurltype)
    case(4)
       curlB(:,i) = rho(i)*curlB(:,i)
    case(3)
       curlB(:,i) = weight*curlB(:,i) !!*gradh(i)
    case(2)
       curlB(:,i) = -rho(i)*curlB(:,i)    
    case default
       curlB(:,i) = curlB(:,i)*gradh(i)/rho(i)
       if (present(curlBgradh)) then
          curlBgradh(:,i) = curlBgradh(:,i)*gradh(i)
       endif
    end select
 enddo

 return
end subroutine get_curl
      
end module getcurl
