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

module density_summations
  use get_neighbour_lists
  implicit none
    
contains
!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab (h_a)
!!
!! Also computes the variable smoothing length terms sum_b m_b dW_ab/dh_a
!!
!! This version computes the density on all particles
!! and therefore only does each pairwise interaction once
!!------------------------------------------------------------------------

  subroutine density(x,pmass,hh,vel,rho,drhodt,densn,dndt, &
                     gradh,gradhn,gradsoft,gradgradh,npart,ntotal)
    use dimen_mhd,    only:ndim,ndimV
    use debug,        only:trace
    use loguns,       only:iprint
    use kernels,      only:radkern2,interpolate_kernel,interpolate_kernels_dens,interpolate_kernel_soft
    use linklist,     only:ll,ifirstincell,numneigh,ncellsloop
    use options,      only:ikernav,igravity,imhd,ikernel,ikernelalt,iprterm
    use part,         only:Bfield,uu,psi,itype,itypebnd
    use setup_params, only:hfact
    use rates,        only:dBevoldt
!
!--define local variables
!
    implicit none
    integer, intent(in) :: npart,ntotal
    real, dimension(:,:), intent(in) :: x, vel
    real, dimension(:), intent(in) :: pmass, hh
    real, dimension(:), intent(out) :: rho,drhodt,densn,dndt,gradh,gradhn,gradsoft,gradgradh
 
    integer :: i,j,n
    integer :: icell,iprev,nneigh
    integer, dimension(ntotal) :: listneigh
    integer :: idone
    integer, parameter :: itemp = 121
!
!  (particle properties - local copies)
!      
    real :: rij,rij2
    real :: hi,hi1,hav,hav1,hj,hj1,hi21
    real :: hfacwab,hfacwabi,hfacwabj
    real, dimension(ndim) :: dx,xi
    real, dimension(ndimV) :: veli,dvel
    real, dimension(ntotal) :: rhoin
    real, dimension(ntotal) :: h1,unity
    real :: dvdotr,pmassi,pmassj,projBi,projBj
    real, dimension(ndimV) :: dr
    integer :: itypei,itypej
!
!  (kernel quantities)
!
    real :: q2,q2i,q2j      
    real :: wab,wabi,wabj,weight,wabalti,wabaltj
    real :: grkern,grkerni,grkernj,grkernalti,grkernaltj,grgrkerni,grgrkernj
    real :: dwdhi,dwdhj,dwaltdhi,dwaltdhj,dphidhi,dphidhj,dwdhdhi,dwdhdhj ! grad h terms
    real :: wconst
!
!--allow for tracing flow
!      
    if (trace) write(iprint,*) ' Entering subroutine density'  
!
!--initialise quantities
!
    listneigh = 0
    numneigh = 0
    dwdhi = 0.
    dwdhj = 0.
    dwaltdhi = 0.
    dwaltdhj = 0.
    dwdhdhi = 0.
    dwdhdhj = 0.
    wconst = 1./hfact**ndim
    dr(:) = 0.

    do i=1,npart
       if (itype(i) /= itypebnd) then
          rhoin(i) = rho(i)
          rho(i) = 0.
          drhodt(i) = 0.
          densn(i) = 0.
          dndt(i) = 0.
          gradh(i) = 0.
          gradhn(i) = 0.
          gradsoft(i) = 0.
          gradgradh(i) = 0.
          if (imhd.eq.5) dBevoldt(:,i) = 0.
          if (imhd.eq.0) then
             psi(i) = 0.
             unity(i) = 0.
          endif
       endif
    enddo
    do i=1,ntotal
       h1(i) = 1./hh(i)
    enddo
!
!--Loop over all the link-list cells
!
    loop_over_cells: do icell=1,ncellsloop   ! step through all cells
!
!--get the list of neighbours for this cell 
!  (common to all particles in the cell)
!
       call get_neighbour_list(icell,listneigh,nneigh)
!
!--now loop over all particles in the current cell
!
       i = ifirstincell(icell)
       idone = -1       ! note density summation includes current particle
       if (i.NE.-1) iprev = i

       loop_over_cell_particles: do while (i.NE.-1)              ! loop over home cell particles

!       PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
          idone = idone + 1
          pmassi = pmass(i)
          xi(:) = x(:,i)
          veli(:) = vel(:,i) 
          hi = hh(i)
          hi1 = h1(i)
          hfacwabi = hi1**ndim
          hi21 = hi1*hi1
          itypei = itype(i)
!
!--for each particle in the current cell, loop over its neighbours
!
          loop_over_neighbours: do n = idone+1,nneigh
             j = listneigh(n)
             !--skip particles of different type
             itypej = itype(j)
             if (itypej.ne.itypei .and. itypej.ne.itypebnd .and. itypei.ne.itypebnd) cycle loop_over_neighbours

             dx(:) = xi(:) - x(:,j)

             hj = hh(j)
             hj1 = h1(j)

             rij2 = dot_product(dx,dx)
             q2i = rij2*hi21
             q2j = rij2*hj1*hj1
!          PRINT*,' neighbour,r/h,dx,hi,hj ',j,SQRT(q2i),dx,hi,hj
!       
!--do interaction if r/h < compact support size
!  don't calculate interactions between ghost particles
!
             if ((q2i.LT.radkern2).OR.(q2j.LT.radkern2)  &
                  .AND. (i.LE.npart .OR. j.LE.npart)) then

                !if (i.eq.itemp .or. j.eq.itemp) then
                !   print*,' neighbour,r/hi,r/hj,hi,hj:',i,j,sqrt(q2i),sqrt(q2j),hi,hj,rho(itemp)
                !endif

                if (i.LE.npart) numneigh(i) = numneigh(i) + 1
                if (j.LE.npart .and. j.ne.i) numneigh(j) = numneigh(j) + 1
                rij = sqrt(rij2)
                dr(1:ndim) = dx(1:ndim)/(rij + epsilon(rij))
                hfacwabj = hj1**ndim         
!
!--weight self contribution by 1/2
!
                if (j.EQ.i) then
                   weight = 0.5
                else
                   weight = 1.0
                endif
                pmassj = pmass(j)
!       
!--interpolate from kernel table              
!  (use either average h or average kernel gradient)
!
                if (ikernav.EQ.1) then              
!  (using average h)
                   hav = 0.5*(hi + hj)
                   hav1 = 1./hav
                   hfacwab = hav1**ndim
                   q2 = rij2*hav1*hav1
                   call interpolate_kernel(q2,wab,grkern)
                   wab = wab*hfacwab
                   grkern = grkern*hfacwab*hav1
                   wabi = wab
                   wabj = wab
                   grkerni = grkern
                   grkernj = grkern
                else
                   if (igravity.ne.0 .and. ikernav.eq.3 .and. ndim.eq.3) then
                      call interpolate_kernel_soft(q2i,wabi,grkerni,dphidhi)
                      gradsoft(i) = gradsoft(i) + weight*pmassj*dphidhi*hi21
                      call interpolate_kernel_soft(q2j,wabj,grkernj,dphidhj)
                      gradsoft(j) = gradsoft(j) + weight*pmassi*dphidhj*hj1*hj1
                      if (ikernelalt.ne.ikernel) then
                         call interpolate_kernels_dens(q2i,wabi,grkerni,grgrkerni,wabalti,grkernalti)
                         call interpolate_kernels_dens(q2j,wabj,grkernj,grgrkernj,wabaltj,grkernaltj)
                      else
                         stop 'grgrkerni will not work here'
                         wabalti = wabi
                         grkernalti = grkerni
                         wabaltj = wabj
                         grkernaltj = grkernj
                         grgrkerni = 0.
                         grgrkernj = 0.
                      endif
                   else
                      call interpolate_kernels_dens(q2i,wabi,grkerni,grgrkerni,wabalti,grkernalti)
                      call interpolate_kernels_dens(q2j,wabj,grkernj,grgrkernj,wabaltj,grkernaltj)
                   endif
              !  (using hi)
                   wabi = wabi*hfacwabi
                   wabalti = wabalti*hfacwabi
                   grkerni = grkerni*hfacwabi*hi1
                   grgrkerni = grgrkerni*hfacwabi*hi1*hi1
                   grkernalti = grkernalti*hfacwabi*hi1
              !  (using hj)
                   wabj = wabj*hfacwabj
                   wabaltj = wabaltj*hfacwabj
                   grkernj = grkernj*hfacwabj*hj1
                   grgrkernj = grgrkernj*hfacwabj*hj1*hj1
                   grkernaltj = grkernaltj*hfacwabj*hj1
              !  (calculate average)
                   if (ikernav.eq.2) then
                      wab = 0.5*(wabi + wabj)
                      wabi = wab
                      wabj = wab
                      grkern = 0.5*(grkerni + grkernj)
                      grkerni = grkern
                      grkernj = grkern
                   endif
              !
              !--derivative w.r.t. h for grad h correction terms (and dhdrho)
              !              
                   dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
                   dwdhj = -rij*grkernj*hj1 - ndim*wabj*hj1
                   dwaltdhi = -rij*grkernalti*hi1 - ndim*wabalti*hi1
                   dwaltdhj = -rij*grkernaltj*hj1 - ndim*wabaltj*hj1
                   
                   dwdhdhi = ndim*(ndim+1)*wabi*hi1**2 + 2.*(ndim+1)*rij*hi1**2*grkerni &
                           + rij**2*hi1**2*grgrkerni
                   dwdhdhj = ndim*(ndim+1)*wabj*hj1**2 + 2.*(ndim+1)*rij*hj1**2*grkernj &
                           + rij**2*hj1**2*grgrkernj
                endif
!
!--calculate density and number density
!
                if (itypei /= itypebnd) then
                   rho(i) = rho(i) + pmassj*wabi*weight
                   densn(i) = densn(i) + wabalti*weight
                endif
                
                if (itypej /= itypebnd) then
                   rho(j) = rho(j) + pmassi*wabj*weight
                   densn(j) = densn(j) + wabaltj*weight
                endif
!
!--drhodt, dndt
!
                if (i.ne.j) then
                   dvel(1:ndimV) = veli(1:ndimV) - vel(1:ndimV,j)
                   dvdotr = dot_product(dvel,dr)
                   drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni !+ pmassj*dvdotr*wabi*hi1
                   drhodt(j) = drhodt(j) + pmassi*dvdotr*grkernj !+ pmassi*dvdotr*wabj*hj1
                   dndt(i) = dndt(i) + dvdotr*grkernalti
                   dndt(j) = dndt(j) + dvdotr*grkernaltj
                   if (imhd.eq.5) then
                   projBi = dot_product(Bfield(:,i),dr)
                   projBj = dot_product(Bfield(:,j),dr)
                   dBevoldt(:,i) = dBevoldt(:,i) - pmassj*projBi*dvel(:)*grkerni
                   dBevoldt(:,j) = dBevoldt(:,j) - pmassi*projBj*dvel(:)*grkernj
                   endif
                else
                   !drhodt(i) = drhodt(i) + pmassj*dvdotr*wabi*hi1
                endif

                if (ikernav.EQ.3) then
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known
!  also do the number density version

                   if (itypei /= itypebnd) then
                      gradh(i) = gradh(i) + weight*pmassj*dwdhi
                      gradhn(i) = gradhn(i) + weight*dwaltdhi
                      gradgradh(i) = gradgradh(i) + weight*pmassj*dwdhdhi
                   endif
                   if (itypej /= itypebnd) then
                      gradh(j) = gradh(j) + weight*pmassi*dwdhj
                      gradhn(j) = gradhn(j) + weight*dwaltdhj
                      gradgradh(j) = gradgradh(j) + weight*pmassi*dwdhdhj
                   endif
                endif
                
                if (imhd.eq.0) then
                   if (iprterm.eq.10) then
                      psi(i) = psi(i) + pmassj*wabi*uu(j)
                      unity(i) = unity(i) + wconst*wabi/hfacwabi
                      if (i.ne.j) then
                         psi(j) = psi(j) + pmassi*wabj*uu(i)
                         unity(j) = unity(j) + wconst*wabj/hfacwabj
                      endif
                   elseif (iprterm.eq.12) then
                      psi(i) = psi(i) + wconst*wabi/hfacwabi
                      if (i.ne.j) then
                         psi(j) = psi(j) + wconst*wabj/hfacwabj
                      endif                   
                   endif
                endif
                
             endif
           
          enddo loop_over_neighbours

          iprev = i
          if (iprev.NE.-1) i = ll(i)  ! possibly should be only IF (iprev.NE.-1)
       enddo loop_over_cell_particles

    enddo loop_over_cells

    if (imhd.eq.5) then
       do i=1,npart
          dBevoldt(:,i) = dBevoldt(:,i)/rhoin(i)**2
       enddo
    endif
    return
  end subroutine density
      

!!------------------------------------------------------------------------
!! Computes the density by direct summation over the particles neighbours
!! ie. rho_a = sum_b m_b W_ab (h_a)
!!
!! This version computes the density only on a selected list of particles
!! given by the contents of the array ipartlist (enables iteration on 
!! unconverged particles only). It is therefore slightly slower
!! since some particle pairs may be done twice, once for each particle.
!!
!! Assumes rho_a is only a function of h_a
!!
!! This version must be used for individual particle timesteps
!!------------------------------------------------------------------------
  
  subroutine density_partial(x,pmass,hh,vel,rho,drhodt,densn,dndt, &
                             gradh,gradhn,gradsoft,gradgradh,ntotal,nlist,ipartlist)
    use dimen_mhd,  only:ndim,ndimV
    use debug,      only:trace
    use loguns,     only:iprint
 
    use kernels,      only:radkern2,interpolate_kernels_dens,interpolate_kernel_soft
    use linklist,     only:iamincell,numneigh
    use options,      only:igravity,imhd,ikernel,ikernelalt,iprterm
    use part,         only:Bfield,uu,psi,itype,itypebnd
    use rates,        only:dBevoldt
    use setup_params, only:hfact
!
!--define local variables
!
    implicit none
    integer, intent(in) :: ntotal
    real, dimension(:,:), intent(in) :: x, vel
    real, dimension(:), intent(in) :: pmass, hh
    real, dimension(:), intent(out) :: rho,drhodt,densn,dndt,gradh,gradhn,gradsoft,gradgradh
    integer, intent(in) :: nlist
    integer, intent(in), dimension(:) :: ipartlist

    integer :: i,j,n
    integer :: icell,ipart,nneigh !!,minneigh,minpart
    integer, dimension(ntotal) :: listneigh
    integer :: icellprev
!
!  (particle properties - local copies)
!      
    real :: rij,rij2
    real :: hi,hi1,hi2,hi21
    real :: hfacwabi,hfacgrkerni,pmassj
    real, dimension(ndim)  :: dx,xi
    real, dimension(ndimV) :: veli,dvel
    real, dimension(nlist) :: rhoin
    real, dimension(ndimV) :: dr
    real :: dvdotr,projBi
    integer :: itypei
!
!  (kernel quantities)
!
    real :: q2i      
    real :: wabi,wabalti,grkerni,grgrkerni,grkernalti 
    real :: dwdhi,dwaltdhi,dphidhi,dwdhdhi ! grad h terms
    real :: wconst,unityi
!
!--allow for tracing flow
!      
    if (trace) write(iprint,*) ' Entering subroutine density_partial'  
!
!--initialise quantities
!
    listneigh = 0
    wconst = 1./hfact**ndim
    dr(:) = 0.

    do ipart=1,nlist
       i = ipartlist(ipart)
       rhoin(ipart) = rho(i)
       rho(i) = 0.
       drhodt(i) = 0.
       densn(i) = 0.
       dndt(i) = 0.
       gradh(i) = 0.
       gradhn(i) = 0.
       gradsoft(i) = 0.
       gradgradh(i) = 0.
       numneigh(i) = 0
       if (imhd.eq.5) dBevoldt(:,i) = 0.
       if (imhd.eq.0 .and. iprterm.eq.10) psi(i) = 0.
    enddo
    icellprev = 0
!
!--Loop over all the particles in the density list
!
    loop_over_particles: do ipart=1,nlist            ! step through all cells

       i = ipartlist(ipart)
!
!--find cell of current particle
!
       icell = iamincell(i)
!    PRINT*,' particle ',i,' cell = ',icell
!
!--if different to previous cell used, get the list of neighbours for this cell
!  (common to all particles in the cell)
!
       if (icell.NE.icellprev) then
          call get_neighbour_list_partial(icell,listneigh,nneigh)
       endif
       icellprev = icell

   !!!PRINT*,'Doing particle ',i,nneigh,' neighbours',hh(i)
       hi = hh(i)
       hi2 = hi*hi
       hi1 = 1./hi
       hi21 = hi1*hi1
       unityi = 0.
       hfacwabi = hi1**ndim
       hfacgrkerni = hfacwabi*hi1
       xi = x(:,i)
       veli(:) = vel(:,i) 
       itypei = itype(i)
       if (itypei.eq.itypebnd) cycle loop_over_particles
!
!--loop over current particle's neighbours
!
       loop_over_neighbours: do n = 1,nneigh
          j = listneigh(n)
          !--skip particles of different type
          if (itype(j).ne.itypei .and. itype(j).ne.itypebnd) cycle loop_over_neighbours

          dx(:) = xi(:) - x(:,j)
!
!--calculate averages of smoothing length if using this averaging
!                           
          rij2 = dot_product(dx,dx)
          q2i = rij2*hi21
!      
!--do interaction if r/h < compact support size
!
          if (q2i.LT.radkern2) then
             rij = sqrt(rij2)
             dr(1:ndim) = dx(1:ndim)/(rij + epsilon(rij))
             !!!if (i.eq.416) PRINT*,' neighbour,r/h,hi ',j,SQRT(q2i),hi
             numneigh(i) = numneigh(i) + 1
             pmassj = pmass(j)
!      
!--interpolate from kernel table (using hi)
!
             if (igravity.ne.0) then
                call interpolate_kernel_soft(q2i,wabi,grkerni,dphidhi)
                gradsoft(i) = gradsoft(i) + pmassj*dphidhi/hi2
                if (ikernel.ne.ikernelalt) then
                   call interpolate_kernels_dens(q2i,wabi,grkerni,grgrkerni,wabalti,grkernalti)             
                else
                   wabalti = wabi
                   grkernalti = grkerni
                endif
             else
                call interpolate_kernels_dens(q2i,wabi,grkerni,grgrkerni,wabalti,grkernalti)             
             endif
             wabi = wabi*hfacwabi
             wabalti = wabalti*hfacwabi
             grkerni = grkerni*hfacgrkerni
             grgrkerni = grgrkerni*hfacwabi*hi1*hi1
             grkernalti = grkernalti*hfacgrkerni
!
!--derivative w.r.t. h for grad h correction terms (and dhdrho)
!
             dwdhi = -rij*grkerni*hi1 - ndim*wabi*hi1
             dwaltdhi = -rij*grkernalti*hi1 - ndim*wabalti*hi1

             dwdhdhi = ndim*(ndim+1)*wabi*hi1**2 + 2.*(ndim+1)*rij*hi1**2*grkerni &
                     + rij**2*hi1**2*grgrkerni

!
!--calculate density and number density
!
             rho(i) = rho(i) + pmassj*wabi
             densn(i) = densn(i) + wabalti
!
!--drhodt, dndt
!
             if (i.ne.j) then
                dvel(:) = veli(:) - vel(:,j)
                dvdotr = dot_product(dvel,dr)
                drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni !+ pmassj*dvdotr*wabi*hi1
                dndt(i) = dndt(i) + dvdotr*grkernalti
                if (imhd.eq.5) then
                projBi = dot_product(Bfield(:,i),dr)
                dBevoldt(:,i) = dBevoldt(:,i) - pmassj*projBi*dvel(:)*grkerni
                endif
             else
                !drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni + pmassj*dvdotr*wabi*hi1
             endif
!
!--correction term for variable smoothing lengths
!  this is the small bit that should be 1-gradh
!  need to divide by rho once rho is known

             gradh(i) = gradh(i) + pmassj*dwdhi
             gradhn(i) = gradhn(i) + dwaltdhi
             gradgradh(i) = gradgradh(i) + pmassj*dwdhdhi

             if (imhd.eq.0 .and. iprterm.eq.10) then
                psi(i) = psi(i) + pmassj*wabi*uu(j)
                unityi = unityi + wconst*wabi/hfacwabi
             endif
          endif
          
       enddo loop_over_neighbours

       if (imhd.eq.0) then
          !psi(i) = abs(uu(i) - psi(i)/unityi)
          !if (psi(i)/uu(i).lt.0.01) psi(i) = 0.
       else
          psi(i) = 0.
       endif

    enddo loop_over_particles

    !print*,'finished density_partial, rho, gradh, h =',rho(1),gradh(1),hh(1),numneigh(1)
    !print*,'maximum number of neighbours = ',MAXVAL(numneigh),MAXLOC(numneigh),rho(MAXLOC(numneigh))
    !print*,'minimum number of neighbours = ',MINVAL(numneigh(1:npart)), &
    !       MINLOC(numneigh(1:npart)),rho(MINLOC(numneigh(1:npart)))

    !minneigh = 100000
    !minpart = 1
    !do i=1,nlist
    !   j = ipartlist(i)
    !   if (numneigh(j).lt.minneigh) then
    !      minneigh = numneigh(j)
    !      minpart = j
    !   endif
    !enddo
    !print*,'minimum number of neighbours = ',minneigh,minpart,rho(minpart)
    if (imhd.eq.5) then
    do i=1,nlist
       j = ipartlist(i)
       dBevoldt(:,j) = dBevoldt(:,j)/rhoin(i)**2
    enddo
    endif

    return
  end subroutine density_partial
      
  
end module density_summations
