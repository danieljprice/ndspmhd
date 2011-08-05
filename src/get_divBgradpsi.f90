!!------------------------------------------------------------------------
!! Computes an SPH estimate of div B and grad psi given B and psi
!! for superfast cleaning
!!
!! This version computes div B on all particles
!! and therefore only does each pairwise interaction once
!!
!! particle quantities are explicitly passed rather than through the 
!! global modules
!!------------------------------------------------------------------------

subroutine get_divBgradpsi(divB,gradpsi,Bin,psi,x,hh,pmass,rho,npart,ntot)
 use dimen_mhd, only:ndim, ndimV
 use debug, only:trace
 use loguns, only:iprint
 
 use bound, only:ireal
 use kernels, only:interpolate_kernel,radkern2
 use linklist
 use options, only:ikernav,imhd
 use get_neighbour_lists
 implicit none
 integer, intent(in) :: npart,ntot
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: hh,pmass
 real, dimension(ndimv,ntot), intent(in) :: bin
 real, dimension(ntot), intent(in) :: psi, rho
 real, dimension(ntot), intent(out) :: divb
 real, dimension(ndimv,ntot), intent(out) :: gradpsi
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
 real :: hfacwab,hfacwabi,hfacwabj,rho1j,rho1i
 real :: pmassi,pmassj,projdb,gradpsiterm
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
 if (trace) write(iprint,*) ' entering subroutine get_divbgradpsi, npart,ntot=',npart,ntot 
!
!--initialise quantities
!
 listneigh = 0
 divb = 0.
 gradpsi = 0.
 !!rho = 0.
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
       rho1i = 1./rho(i)
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
             rho1j = 1./rho(j)
             
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
                gradpsiterm = (psi(i)-psi(j))*grkern ! (-ve grad psi)

                divb(i) = divb(i) - pmassj*projdb*grkerni
                divb(j) = divb(j) - pmassi*projdb*grkernj            
                gradpsi(:,i) = gradpsi(:,i) + pmassj*gradpsiterm*dr(:)
                gradpsi(:,j) = gradpsi(:,j) + pmassi*gradpsiterm*dr(:)
                   
!                vsig = 0.5*alphasub*(spsound(i) + spsound(j) &
!                          + dot_product(bin(:,i),bin(:,i))*rho1i &
!                          + dot_product(bin(:,j),bin(:,j))*rho1j)
!                bdiss(:) = (bin(:,i) - bin(:,j)) - 3.*projdb*dr(:)
!                gradpsi(:,i) = gradpsi(:,i) - vsig*pmassj*rho1j*bdiss(:)*grkern
!                gradpsi(:,j) = gradpsi(:,j) + vsig*pmassi*rho1i*bdiss(:)*grkern
                !print*,' bi,j = ',bin(:,i),bin(:,j)
                !print*,' projdb,dr = ',projdb,dr(:)
                !read*
                   
                !weight = 1.0
                !f (j.eq.i) weight = 0.5
                !rho(i) = rho(i) + weight*pmassj*wabi
                !rho(j) = rho(j) + weight*pmassi*wabj

                !      else
                !         print*,' r/h > 2 '      
                
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
       if (imhd.ge.11) then ! evolving b
          gradpsi(:,i) = gradpsi(:,i)/rho(i)
       else                 ! evolving b/rho
          gradpsi(:,i) = gradpsi(:,i)/rho(i)**2       
       endif
       divb(i) = divb(i)/rho(i)
    else
       write(*,*) 'error: get_divbgradpsi: rho < 0.'
    endif
    !print*,'divb, rho = ',divb(i),rho(i),bin(:,i)
    !if (mod(i,20).eq.0) read*
 enddo
 
 do i=npart+1,ntot
    j = ireal(i)
    !!rho(i) = rho(j)
    gradpsi(:,i) = 0.
    !print*,i,bin(:,i),rho(i)
 enddo
! read*

 return
end subroutine get_divbgradpsi
      
