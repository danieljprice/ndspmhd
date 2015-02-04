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

module dustdiffusion
 implicit none
 
contains
!!------------------------------------------------------------------------
!! Computes dust diffusion term for one fluid dust in the terminal
!! velocity approximation (assumes deltav is known)
!!
!! particle quantities are explicitly passed rather than through the 
!! global modules
!!
!! assumes there has been a previous call to density
!!------------------------------------------------------------------------

subroutine dust_diffusion(npart,ntot,x,pmass,rho,hh,gradh,dustfrac,ddustevoldt,deltav,vel,pr,uu,spsound,dudt)
 use dimen_mhd, only:ndim, ndimV
 use debug,     only:trace
 use loguns,    only:iprint
 use bound,     only:ireal
 use kernels,   only:interpolate_kernel,radkern2
 use linklist
 use options,   only:iener,Kdrag,use_sqrtdustfrac,idrag_nature
 use dust,      only:get_tstop
 use get_neighbour_lists
 integer, intent(in) :: npart,ntot
 real, dimension(ndim,ntot),  intent(in) :: x
 real, dimension(ntot),       intent(in) :: pmass,rho,hh,gradh,dustfrac,uu,spsound,pr
 real, dimension(ndimV,ntot), intent(in) :: vel
 real, dimension(ndimV,ntot), intent(inout) :: deltav
 real, dimension(ntot),      intent(out) :: ddustevoldt
 real, dimension(ntot),    intent(inout) :: dudt
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
 real :: hi,hi1,hj,hj1
 real :: hfacwabi,hfacwabj
 real :: pmassi,pmassj,rhoi,rhoj,gradhi,rho1i,rho1j,rhogasi,rhogasj
 real :: dustfraci,dustfracj,projdeltavi,projdeltavj,dvgasdotr
 real :: termi,termj,term,du,disstermi,disstermj,deltav2i,deltav2j
 real :: tstopi,tstopj,pdvtermi,pdvtermj,si,sj,rhodusti,rhodustj
 real, dimension(ndim) :: dx
 real, dimension(ndimv) :: dr,deltavi,deltavj,dvgas,vgasi,vgasj
!
!  (kernel quantities)
!
 real :: q2i,q2j
 real :: wabi,wabj
 real :: grkerni,grkernj
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine dust_diffusion'
!
!--initialise quantities
!  (do NOT set du/dt to zero)
!
 !ddustevoldt = 0.
 listneigh = 0
!
! make sure deltav has been copied to ghosts
!
 do i=npart+1,ntot
    j = ireal(i)
    deltav(:,i) = deltav(:,j)
 enddo
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

       idone = idone + 1
       hi = hh(i)
       hi1 = 1./hi
       hfacwabi = hi1**ndim
       pmassi    = pmass(i)
       rhoi      = rho(i)
       rho1i     = 1./rhoi
       gradhi    = gradh(i)
       dustfraci = dustfrac(i)
       rhogasi   = (1. - dustfraci)*rhoi
       rhodusti  = rhoi - rhogasi
       deltavi   = deltav(:,i)
       deltav2i  = dot_product(deltavi,deltavi)
       vgasi     = vel(:,i) !- dustfraci*deltavi
       tstopi    = get_tstop(idrag_nature,rhogasi,rhodusti,spsound(i),Kdrag)
       !tstopi    = dustfraci*rhogasi/Kdrag
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
          j = listneigh(n)
          if ((j.ne.i) .and..not.(j.gt.npart .and. i.gt.npart)) then
             ! don't count particle with itself       
             dx(:) = x(:,i) - x(:,j)
             hj = hh(j)
             hj1 = 1./hj
             hfacwabj = hj1**ndim

             rij2 = dot_product(dx,dx)             
             q2i = rij2*hi1*hi1
             q2j = rij2*hj1*hj1
!     
!--do interaction if r/h < compact support size
!
             if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then
                rij = sqrt(rij2)
                pmassj    = pmass(j)
                rhoj      = rho(j)
                rho1j     = 1./rhoj
                dustfracj = dustfrac(j)
                rhogasj   = (1. - dustfracj)*rhoj
                rhodustj  = rhoj - rhogasj
                deltavj   = deltav(:,j)
                deltav2j  = dot_product(deltavj,deltavj)
                tstopj    = get_tstop(idrag_nature,rhogasj,rhodustj,spsound(j),Kdrag)
                !tstopj    = dustfracj*rhogasj/Kdrag
                dr = 0.
                dr(1:ndim) = dx(1:ndim)/(rij + epsilon(rij))  ! unit vector
!     
!--interpolate from kernel table
!
                !  (using hi)
                call interpolate_kernel(q2i,wabi,grkerni)
                wabi = wabi*hfacwabi
                grkerni = grkerni*hfacwabi*hi1*gradhi
                !  (using hj)
                call interpolate_kernel(q2j,wabj,grkernj)
                wabj = wabj*hfacwabj
                grkernj = grkernj*hfacwabj*hj1*gradh(j)
!
!--calculate depsilon/dt
!
                projdeltavi = dot_product(deltav(:,i),dr)
                projdeltavj = dot_product(deltav(:,j),dr)

                if (use_sqrtdustfrac) then
                   si = sqrt(rhoi*dustfraci)
                   sj = sqrt(rhoj*dustfracj)
                   termi = (1. - dustfraci)*rho1i*projdeltavi*grkerni
                   termj = (1. - dustfracj)*rho1j*projdeltavj*grkernj
                else
                   termi = dustfraci*(1. - dustfraci)*rho1i*projdeltavi*grkerni
                   termj = dustfracj*(1. - dustfracj)*rho1j*projdeltavj*grkernj
                   si = 1.
                   sj = 1.
                endif
                term = termi + termj
                ddustevoldt(i) = ddustevoldt(i) - pmassj*sj*term
                ddustevoldt(j) = ddustevoldt(j) + pmassi*si*term

                if (iener.gt.0) then
                   disstermi = 0.!dustfraci*deltav2i/tstopi
                   disstermj = 0.!dustfracj*deltav2j/tstopj
                   
                   vgasj = vel(:,j) !- dustfracj*deltavj
                   dvgas = vgasi - vgasj
                   dvgasdotr = dot_product(dvgas,dr)

                   pdvtermi = pr(i)*rho1i/rhogasi*pmassj*dvgasdotr*grkerni
                   pdvtermj = pr(j)*rho1j/rhogasj*pmassi*dvgasdotr*grkernj

                   if (iener.ne.2) stop 'only thermal energy equation implemented for idust=3'
                   du = uu(i) - uu(j)
!                   dudt(i) = dudt(i) - pmassj*termi/(1. - dustfraci)*du
!                   dudt(j) = dudt(j) - pmassi*termj/(1. - dustfracj)*du
                   dudt(i) = dudt(i) - dustfraci*rho1i*pmassj*du*projdeltavi*grkerni + pdvtermi + disstermi
                   dudt(j) = dudt(j) - dustfracj*rho1j*pmassi*du*projdeltavj*grkernj + pdvtermj + disstermj
                endif
             endif

         endif! .not. j>npart .and. i>npart
       enddo loop_over_neighbours

       iprev = i
       if (iprev.ne.-1) i = ll(i)
    enddo loop_over_cell_particles
    
 enddo loop_over_cells

 do i=npart+1,ntot
    j = ireal(i)
    ddustevoldt(i) = ddustevoldt(j)
    dudt(i) = dudt(j)
 enddo
 
 return
end subroutine dust_diffusion
      
end module dustdiffusion
