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

!!--------------------------------------------------------------------
!! Computes the rates of change of the conserved variables
!! (forces, energy etc)
!! This is the core of the SPH algorithm
!!--------------------------------------------------------------------

subroutine get_rates
! USE dimen_mhd
 use debug,  only:trace
 use loguns, only:iprint
 use bound,  only:pext
 use artvi
 use eos
 use hterms
 use kernels, only:radkern2
 use linklist
 use options
 use part
 use rates
 use timestep
 use xsph

 use fmagarray
 use derivB
 use get_neighbour_lists
 use grutils,       only:metric_diag,dot_product_gr
 use getBeulerpots, only:compute_rmatrix,exactlinear
 use matrixcorr,    only:dxdx,idxdx,jdxdx,ndxdx
 use utils,         only:cross_product3D
 use resistivity,   only:etafunc
 use externf,       only:external_forces,pequil
 use dust,          only:get_tstop
!
!--define local variables
!
 implicit none
 integer :: i,j,n,k
 integer :: icell,iprev,nneigh,nlistdim
 integer, allocatable, dimension(:) :: listneigh
 integer :: idone,nclumped
!
!  (particle properties - local copies and composites)
!
 real :: rij,rij2
 real :: rhoi,rho1i,rho2i,rho21i,rhoj,rho1j,rho2j,rho21j,rhoav1,rhoij
 real :: pmassi,pmassj,projvi,projvj
 real :: Prho2i,Prho2j,prterm,pri,prj
 real :: hi,hi1,hj,hj1,hi21,hj21
 real :: hfacwabi,hfacgrkerni
 real, dimension(ndim) :: xi, dx
 real, dimension(ndimV) :: fexternal  !!,dri,drj
 real, dimension(ntotal) :: h1
 integer :: itypei,itypej
!
!--gr terms
!
 real :: sqrtgi,sqrtgj,dens1i,v2i,v2j
 real, dimension(ndimV) :: gdiagi,gdiagj
!
!  (velocity)
!
 real, dimension(ndimV) :: veli,velj,dvel,fmean
 real, dimension(ndimV) :: dr
 real :: dvdotr
!
!  (mhd)
!
 real, dimension(ndimB) :: Brhoi,Brhoj,Bi,Bj,dB
 real, dimension(ndimB) :: faniso,fmagi,fmagj,Bevoli,dBevol
 real, dimension(ndimB) :: curlBi,forcei,forcej,dBevoldti
 real :: fiso,B2i,B2j
 real :: valfven2i,valfven2j
 real :: BidotdB,BjdotdB,Brho2i,Brho2j,BdotBexti
 real :: projBrhoi,projBrhoj,projBi,projBj,projdB,projBconst
 real, dimension(:,:), allocatable :: curlBsym
 real, dimension(:), allocatable :: divBsym
 real, dimension(:,:,:), allocatable :: dveldx
 real, dimension(6) :: rmatrix
 real :: denom,ddenom,etai,etaj
!
!  (artificial viscosity quantities)
!
 real :: vsig,vsigi,vsigj,vsigav
 real :: spsoundi,spsoundj,alphai,alphaui,alphaBi
!! real :: rhoi5,rhoj5
 real :: vsig2i,vsig2j,vsigproji,vsigprojj,vsigB,vsigu
!! real :: vsigii,vsigjj
 real :: prneti,prnetj
!
!  (av switch)
!
 real :: source,tdecay1,sourcedivB,sourceJ,sourceB,sourceu
 real :: graduterm, graddivvmag,curr2
 real, dimension(:), allocatable :: del2u
!
!  (alternative forms)
!
 real, dimension(:), allocatable :: phi
 real :: phii,phii1,phii_on_phij,phij_on_phii,uui,uuj
!
!  (one fluid dust)
!
 real, dimension(ndimV) :: deltavi,deltavj
 real :: dustfraci(ndust),dustfracj(ndust)
 real :: rhogrhodonrhoi(ndust), rhogrhodonrhoj(ndust)
 real :: deltav2i,deltav2j
 real :: dtstop,projdeltavi,projdeltavj
 real :: rhogasi,rhogasj,projdvgas
 real :: rhodusti(ndust),rhodustj(ndust)
 real, dimension(ndimV) :: vgasi,vgasj,dvgas,fextrai,fextraj
 real :: esum,tstop,ratio,dvmax

!
!  (kernel related quantities)
!
 real :: q2i,q2j
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj
 real :: gradhi,gradhni
 real :: grkernalti,grkernaltj,grgrkernalti,grgrkernaltj
!
!  (time step criteria)
!
 real :: vsigdtc,zero,fhmax, fonh, forcemag, ts_min, h_on_csts_max
 integer :: ierr
!
!--div B correction
!
 real :: gradpsiterm,vsig2,vsigmax !!,dtcourant2
 real :: stressmax,stressterm
 logical, parameter :: itiming = .false.
 real :: t1,t2,t3,t4,t5
 real, allocatable :: fprev(:,:)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine get_rates'

 if (itiming) call cpu_time(t1)
!
!--allocate memory for local arrays
!
 nlistdim = ntotal
 allocate ( listneigh(nlistdim),STAT=ierr )
 if (ierr.ne.0) write(iprint,*) ' Error allocating neighbour list, ierr = ',ierr
 allocate ( phi(ntotal), del2u(ntotal), STAT=ierr )
 if (ierr.ne.0) write(iprint,*) ' Error allocating phi, ierr = ',ierr
 if (imhd.lt.0) then
    allocate( curlBsym(ndimV,ntotal), STAT=ierr )
    if (ierr.ne.0) write(iprint,*) ' Error allocating curlBsym, ierr = ',ierr
    allocate( divBsym(ntotal), STAT=ierr )
    if (ierr.ne.0) write(iprint,*) ' Error allocating divBsym, ierr = ',ierr
 endif
 if (imhd.ne.0 .and. iuse_exact_derivs.gt.0) then
    allocate( dveldx(ndim,ndimV,ntotal), STAT=ierr )
    if (ierr.ne.0) write(iprint,*) ' Error allocating dveldx, ierr = ',ierr
 endif
 listneigh = 0
!
!--initialise quantities
!
 dtcourant  = 1.e6
 dtav       = huge(dtav)
 dtvisc     = huge(dtvisc)
 ts_min     = huge(ts_min)
 h_on_csts_max = 0.
 zero       = 1.e-10
 vsigmax    = 0.
 dr(:)      = 0.
 nclumped   = 0

 do i=1,ntotal      ! using ntotal just makes sure they are zero for ghosts
  force(:,i) = 0.0
  dudt(i) = 0.0
  dendt(i) = 0.0
  dBevoldt(:,i) = 0.0
  daldt(:,i) = 0.0
  dpsidt(i) = 0.0
  gradpsi(:,i) = 0.0
  fmag(:,i) = 0.0
  divB(i) = 0.0
  if (imhd.gt.0 .and. iambipolar.eq.0 .and. iresist.ne.4) curlB(:,i) = 0.0
  if (allocated(curlBsym)) curlBsym(:,i) = 0.
  if (allocated(divBsym)) divBsym(i) = 0.
  if (allocated(dveldx)) then
     dveldx(:,:,i) = 0.
     dxdx(:,i) = 0.
  endif
  if (allocated(del2v)) del2v(:,i) = 0.
  xsphterm(:,i) = 0.0
  del2u(i) = 0.0
  graddivv(:,i) = 0.0
  h1(i) = 1./hh(i)
  if (icty.ge.1) then
     drhodt(i) = 0.
     gradh(i) = 1.
     gradhn(i) = 0.
     dhdt(i) = 0.
  endif
  if (onef_dust) then
     ddustevoldt(:,i) = 0.
     ddeltavdt(:,i) = 0.
  endif
 enddo

!
!--calculate maximum neg stress for instability correction
!
 stressmax = 0.
 if (imhd.ne.0 .and. (imagforce.eq.2 .or. imagforce.eq.7)) then
    do i=1,ntotal
       call metric_diag(x(:,i),gdiagi(:),sqrtgi,ndim,ndimV,geom)
       B2i = dot_product_gr(Bfield(:,i),Bfield(:,i),gdiagi(:))
       if (imagforce.eq.7) then ! vec potential force
          !stressterm = max(1.5*B2i - pr(i),0.)
          !stressmax = max(stressterm,stressmax)
       else
          stressterm = max(0.5*B2i - pr(i),0.)
          stressmax = max(stressterm,stressmax,maxval(Bconst(:)))
       endif
    enddo
    if (stressmax.gt.tiny(stressmax)) write(iprint,*) 'stress correction = ',stressmax
 endif
!
!--set MHD quantities to zero if mhd not set
!
 if (imhd.eq.0) then  ! these quantities are still used if mhd off
    Bi(:) = 0.
    Bj(:) = 0.
    Brhoi(:) = 0.
    Brhoj(:) = 0.
    Brho2i = 0.
    Brho2j = 0.
    valfven2i = 0.
    valfven2j = 0.
    projBi = 0.
    projBj = 0.
    projBrhoi = 0.
    projBrhoj = 0.
    alphaBi = 0.
 endif

!
!--skip the whole neighbour thing if it is doing nothing
!
 if (iprterm.lt.0 .and. iav.eq.0 .and. imhd.eq.0 .and. iener.lt.3 .and. iprterm.ne.-2) then
    write(iprint,*) 'skipping rates'
    goto 666
 endif
!
!  set alternative forms for the SPH equations here
!  phi can be any scalar variable
!
 select case(iprterm)
    case(1)            ! this gives the (P_a + P_b) / (rho_a rho_b) form
       phi(1:ntotal) = rho(1:ntotal)
    case(2)            ! this gives the HK89 form SQRT(Pa Pb)/rhoa rhob
       phi(1:ntotal) = sqrt(pr(1:ntotal))/rho(1:ntotal)
    case(3)
       phi(1:ntotal) = 1./rho(1:ntotal)
    case(4)
       phi(1:ntotal) = 1./rho(1:ntotal)**gamma
    case(5)
       phi(1:ntotal) = rho(1:ntotal)**gamma
    case(6)
       phi(1:ntotal) = 1./en(1:ntotal)
    case(7)
       phi(1:ntotal) = 1./uu(1:ntotal)
    case(8)
       phi(1:ntotal) = 1./pr(1:ntotal)
    case(9)
       phi(1:ntotal) = sqrt(rho(1:ntotal))
    case default      ! this gives the usual continuity, momentum and induction eqns
       phi(1:ntotal) = 1.0
 end select

 if (itiming) call cpu_time(t2)

!
!--Loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop            ! step through all cells

    !!print*,'> doing cell ',icell
    !
    !--get the list of neighbours for this cell
    !  (common to all particles in the cell)
    !
    call get_neighbour_list(icell,listneigh,nneigh)

    i = ifirstincell(icell)      ! start with first particle in cell
    idone = -1
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1)      ! loop over home cell particles
       !!print*,'Doing particle ',i,'of',npart,x(:,i),rho(i),hh(i)
       idone = idone + 1
       xi(:) = x(:,i)
       itypei = itype(i)
       rhoi = rho(i)
       rho2i = rhoi*rhoi
!!       rhoi5 = sqrt(rhoi)
       rho1i = 1./rhoi
       rho21i = rho1i*rho1i
       dens1i = 1./dens(i)
       pri = max(pr(i) - pext,0.)
!       print *, 'PRI ', pr(i), pext, itype(i)
       prneti = pri - pequil(iexternal_force,xi(:),rhoi)
       pmassi = pmass(i)
       Prho2i = pri*rho21i
       spsoundi = spsound(i)
       uui = uu(i)
       veli(:) = vel(:,i)
       v2i = dot_product(veli(:),veli(:))
       alphai = alpha(1,i)
       alphaui = alpha(2,i)
       alphaBi = alpha(3,i)
       phii = phi(i)
       phii1 = 1./phii
       sqrtgi = sqrtg(i)
       ! one fluid dust definitions
       if (onef_dust) then
          dustfraci  = dustfrac(:,i)
          if (use_smoothed_rhodust) then
             rhodusti = rhodust(:,i)
             rhogasi  = rhogas(i)
          else
             rhodusti = dustfraci*rhoi
             rhogasi  = (1. - sum(dustfraci))*rhoi
          endif
          deltavi(:)     = deltav(:,i)
          deltav2i       = dot_product(deltavi,deltavi)
          rhogrhodonrhoi(:) = rhogasi*rhodusti(:)*rho1i
       else
          rhogasi  = rhoi
          rhodusti = 0.
          deltav2i = 0.
       endif
       ! mhd definitions
       if (imhd.ne.0) then
          Bi(:) = Bfield(:,i)
          Brhoi(:) = Bi(:)*rho1i
          BdotBexti = dot_product(Bi(:),Bconst(:))
!          if (geom(1:6).ne.'cartes') then
             call metric_diag(x(:,i),gdiagi(:),sqrtgi,ndim,ndimV,geom)
             B2i = dot_product_gr(Bi,Bi,gdiagi)
!          else
!             B2i = dot_product(Bi,Bi)
!          endif
          Brho2i = B2i*rho21i
          valfven2i = B2i*rho1i
          if (imhd.lt.0) Bevoli(:) = Bevol(:,i)
          if (iresist.eq.3) then
             etai = etafunc(x(1,i),etamhd)
          elseif (iresist > 0) then
             etai = etamhd
          endif
       endif
       gradhi = gradh(i)
       gradhni = gradhn(i)
       hi = hh(i)
       if (hi.le.0.) then
          write(iprint,*) ' rates: h <= 0 particle',i,hi
          call quit
       endif
       hi1 = h1(i)
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       hfacgrkerni = hfacwabi*hi1
       forcei(:) = 0.
       fextrai(:) = 0.
       dBevoldti(:) = 0.
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh

           j = listneigh(n)
           if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then ! don't count particle with itself
           dx(:) = xi(:) - x(:,j)
           !print*,' ... neighbour, h=',j,hh(j),rho(j),x(:,j)
           hj = hh(j)
           hj1 = h1(j) !!1./hj
           hj21 = hj1*hj1

           rij2 = dot_product(dx,dx)
           q2i = rij2*hi21
           q2j = rij2*hj21

             !----------------------------------------------------------------------------
             !  do pairwise interaction if either particle is within range of the other
             !----------------------------------------------------------------------------
           if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then  ! if < 2h
              rij = sqrt(rij2)
              if (rij.le.epsilon(rij) .and. itype(j).eq.itypei) then
                 nclumped = nclumped + 1
                 if (rij.lt.tiny(rij)) then
                    dr = 0.
                    if (itypei /= 2) then
                       write(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj,' type=',itypei
                       call quit
                    endif
                 endif
              elseif (rij.le.epsilon(rij)) then
                 dr = 0.  ! can happen with gas/dust if particles on top of each other
              else
                 dr(1:ndim) = dx(1:ndim)/rij      ! unit vector
              endif
              !do idim=1,ndim
              !   dri(idim) = dot_product(1./gradmatrix(idim,1:ndim,i),dr(1:ndim))
              !   drj(idim) = dot_product(1./gradmatrix(idim,1:ndim,j),dr(1:ndim))
              !enddo
              itypej = itype(j)
              if ((itypej.eq.itypei) .or. (itypei.eq.itypegas .and. itypej.eq.itypebnd) &
                                     .or. (itypej.eq.itypegas .and. itypei.eq.itypebnd) &
                                     .or. (itypei.eq.itypedust .and. itypej.eq.itypebnddust) &
                                     .or. (itypej.eq.itypedust .and. itypei.eq.itypebnddust)) then
!              if (itype(j).eq.itypei &
!                  .or.(itype(j).eq.itypebnd .or. itype(j).eq.itypebnd2) &
!                  .or.(itypei  .eq.itypebnd .or. itype(j).eq.itypebnd2)) then
                 call rates_core
              elseif (idust.eq.2 .and. idrag_nature.gt.0) then !-- drag step if required
                 call drag_forces
              endif
           else      ! if outside 2h
!              PRINT*,'outside 2h, not calculated, r/h=',sqrt(q2i),sqrt(q2j)
           endif      ! r < 2h or not

         endif            ! j.ne.i

        enddo loop_over_neighbours
       !
       !--add contributions to particle i from summation over j
       !  (forcei is forces on gas only, fextra is terms that apply to total fluid)
       !
       force(:,i) = force(:,i) + fextrai(:) + forcei(:)
       dBevoldt(:,i) = dBevoldt(:,i) + dBevoldti(:)
       if (idust.eq.1) ddeltavdt(:,i) = ddeltavdt(:,i) - rhoi/rhogasi*forcei(:)

       iprev = i
       if (iprev.ne.-1) i = ll(i) ! possibly should be only IF (iprev.NE.-1)

    enddo loop_over_cell_particles

 enddo loop_over_cells
 if (nclumped.gt.0) write(iprint,*) ' WARNING: clumping on ',nclumped,' pairs'

666 continue

 if (itiming) call cpu_time(t3)

!----------------------------------------------------------------------------
!  calculate gravitational force on all the particles
!----------------------------------------------------------------------------

 if (igravity.ne.0) then
    fprev = force
    select case(igravity)
    case(1) ! fixed plummer
       call direct_sum_poisson(x(1:ndim,1:npart),pmass(1:npart), &
            potengrav,force(1:ndim,1:npart),hsoft,npart)
    case(2) ! fixed cubic spline
       phi(1:npart) = hsoft
       call direct_sum_poisson_soft(x(1:ndim,1:npart),pmass(1:npart), &
            phi(1:npart),poten(1:npart),force(1:ndim,1:npart),potengrav,npart)
    case(3) ! adaptive cubic spline (no extra term)
       call direct_sum_poisson_soft(x(1:ndim,1:npart),pmass(1:npart), &
            hh(1:npart),poten(1:npart),force(1:ndim,1:npart),potengrav,npart)
    case(4:) ! adaptive cubic spline with energy conservation
       call direct_sum_poisson_soft(x(1:ndim,1:npart),pmass(1:npart), &
            hh(1:npart),poten(1:npart),force(1:ndim,1:npart),potengrav,npart)
    case default
        stop 'unknown igravity setting in forces'
    end select

    !
    !--add v.fgrav term if evolving total energy
    !
    !if (iener==3) then
      ! do i=1,npart
      !    dendt(i) = dendt(i) + dot_product(vel(:,i),force(:,i)-fprev(:,i))
      ! enddo
    !endif
 endif
 if (itiming) call cpu_time(t4)

 if (trace) write(iprint,*) 'Finished main rates loop'

!----------------------------------------------------------------------------
!  loop over the particles again, subtracting external forces and source terms
!----------------------------------------------------------------------------

!
!--calculate maximum vsig over all the particles for use in the hyperbolic cleaning
!
 if (imhd.ne.0 .and. idivBzero.ge.2) then
    vsig2max = vsigmax**2
 endif

 fhmax = 0.0
 fmean(:) = 0.
 dtdrag = huge(dtdrag)
 dtforce= huge(dtforce)
 esum = 0.
 ratio = 0.
 if (idust.eq.2 .and. h_on_csts_max > 1.) then
    print*,' WARNING: violating h < cs*ts resolution criterion by factor of ',h_on_csts_max
 endif

 do i=1,npart

    rhoi  = rho(i)
    rho1i = 1./rhoi
    if (ivisc.gt.0 .and. idust.ne.1 .and. idust.ne.4) then
       esum = esum + pmass(i)*dot_product(vel(:,i),force(:,i))
    endif
    fexternal(:) = 0.
!
!--Dust
!
    if (idust.eq.2 .and. idrag_nature.ne.0 .and. (Kdrag > 0. .or. idrag_nature > 1)) then
       !
       !--two fluid dust: calculate drag timestep
       !
       dtdrag = min(dtdrag,ts_min)
    elseif (idust.eq.1) then
       !------------------
       !  one fluid dust
       !------------------
       if (use_smoothed_rhodust) then
          rhodusti = rhodust(:,i)
          rhogasi  = rhogas(i)
       else
          rhodusti = rhoi*dustfrac(:,i)
          rhogasi  = (1. - sum(dustfrac(:,i)))*rhoi
       endif

       do k=1,ndust
          tstop = get_tstop(idrag_nature,rhogasi,rhodusti(k),spsound(i),Kdrag)
          dtdrag = min(dtdrag,tstop)
          !
          !--d/dt(deltav)  : add terms that do not involve sums over particles
          !
          if (dustfraci(k).gt.0.) then
             dtstop   = 1./tstop
             ddeltavdt(:,i) = ddeltavdt(:,i) - deltav(:,i)*dtstop
          else
             dtstop = 0.
             ddeltavdt(:,i) = 0.
          endif
       enddo

       !
       !--du/dt: add thermal energy dissipation from drag
       !  (maybe should do this in the step routine along with the implicit drag in ddeltav/dt?)
       !
       if (iener.gt.0) then
          deltav2i = dot_product(deltav(:,i),deltav(:,i))
          dudt(i) = dudt(i) + rhodusti(1)*rho1i*deltav2i*dtstop
       endif
    elseif (idust.eq.3 .or. idust.eq.4) then
       !--------------------------------------------
       !  one fluid dust in diffusion approximation
       !--------------------------------------------
       if (use_smoothed_rhodust) then
          rhodusti(:) = rhodust(:,i)
          rhogasi     = rhogas(i)
       else
          rhodusti(:) = rhoi*dustfrac(:,i)
          rhogasi     = (1 - sum(dustfrac(:,i)))*rhoi
       endif
       !
       !--compute stopping time for drag timestep
       !
       dustfraci = rhodusti*rho1i
       do k=1,ndust
          tstop = get_tstop(idrag_nature,rhogasi,rhodusti(k),spsound(i),Kdrag)
          ! CAUTION: Line below must be done BEFORE external forces have been applied
          if (rhogasi > tiny(rhogasi)) then
             deltav(:,i) = -rhoi/rhogasi*force(:,i)*tstop
          else
             !stop 'rhogas < 0'
             deltav(:,i) = 0.
          endif
          ratio = max(dustfraci(k)*tstop/dtcourant,ratio)
          dtdrag = min(dtdrag,0.5*hh(i)**2/(dustfraci(k)*tstop*spsound(i)**2))
       enddo
       if (dtdrag < 0.) then
          !print*,'WARNING: dtdrag = ',dtdrag,dustfraci,dustfrac(:,i),rhodust(:,i),rhogas(i)
          dtdrag = abs(dtdrag)
       endif
       !dvmax = maxval(abs(deltav(:,i)))
       !if (dvmax > 0.) then
       !   dtdrag = min(dtdrag,0.1*hh(i)/dvmax)
       !endif
       if (idrag_nature==4) then
          ! this is constant ts but keeping the particles fixed
          force(:,i) = 0.
          vel(:,i) = 0.
       !   !print*,'dtdrag = ',dtdrag
       endif
    endif

!
!--compute JxB force from the Euler potentials (if using second derivs)
!  also divide by rho for J and div B calculated the "normal" way
!
    if (imhd.ne.0) then
       if (iambipolar > 0) then
          B2i = dot_product(Bfield(:,i),Bfield(:,i))
          dtdrag = min(dtdrag,gamma_ambipolar*rho_ion*rhoi*hh(i)**2/B2i)
       else
          if (imagforce.eq.4) then
             call cross_product3D(curlB(:,i),Bfield(:,i),fmagi(:)) ! J x B
             force(:,i) = force(:,i) + fmagi(:)*rho1i  ! (J x B)/rho
             fmag(:,i) = fmagi(:)*rho1i
          elseif (imhd.gt.0 .and. .not.(iresist==4)) then ! not for vector potential
             curlB(:,i) = curlB(:,i)*rho1i
          endif
       endif
       divB(i) = divB(i)*rho1i

       if (allocated(divBsym)) then
          divBsym(i) = divBsym(i)*rhoi
          divB(i) = divBsym(i)
       endif
    endif
!
!--add external (body) forces
!
    if (iexternal_force.ne.0) then
       call external_forces(iexternal_force,x(1:ndim,i),fexternal(1:ndimV), &
                            ndim,ndimV,vel(1:ndimV,i),hh(i),spsound(i),itype(i))
       force(1:ndimV,i) = force(1:ndimV,i) + fexternal(1:ndimV)
    endif
!
!--add source terms (derivatives of metric) to momentum equation
!
    if (allocated(sourceterms)) force(:,i) = force(:,i) + sourceterms(:,i)
!
!--make dhdt if density is not being done by summation
!  (otherwise this is done in iterate_density)
!
    if (icty.ge.1) then
       dhdt(i) = -hh(i)/(ndim*rho(i))*drhodt(i)
    endif

!
!--calculate maximum force/h for the timestep condition
!  also check for errors in the force
!
!    if ( any(force(:,i).gt.1.e8)) then
!       write(iprint,*) 'rates: force ridiculous ',force(:,i),' particle ',i
!       call quit
!    endif
    fmean(:) = fmean(:) + pmass(i)*force(:,i)
    forcemag = sqrt(dot_product(force(:,i),force(:,i)))
    fonh = forcemag/hh(i)
    if (fonh.gt.fhmax .and. itype(i).ne.1) fhmax = fonh
!
!--calculate simpler estimate of vsig for divergence cleaning and
!  in the dissipation switches
!
    valfven2i = 0.
    if (imhd.ne.0) then
!       if (geom(1:6).ne.'cartes') then
          call metric_diag(x(:,i),gdiagi,sqrtgi,ndim,ndimv,geom)
          valfven2i = dot_product_gr(Bfield(:,i),Bfield(:,i),gdiagi)/dens(i)
!       else
!          valfven2i = dot_product(Bfield(:,i),Bfield(:,i))/dens(i)
!       endif
    endif
    vsig2 = spsound(i)**2 + valfven2i     ! approximate vsig only
    vsig = SQRT(vsig2)
    !!!dtcourant2 = min(dtcourant2,hh(i)/vsig)

    if (iuse_exact_derivs.gt.0) then
       call compute_rmatrix(dxdx(:,i),rmatrix,denom,ndim)
      ! print*,i,' dxdx=  ',dxdx(:,i)
      ! print*,i,' rmatrix = ',rmatrix*rho1i*rho1i,'denom= ',denom*rho1i*rho1i*rho1i
      ! print*,i,' normal derivs, dvxdy = ',dveldx(2,1,i)*rho1i
       if (abs(denom).gt.epsilon(denom)) then
          ddenom = 1./denom
          do k=1,ndimV
             call exactlinear(dveldx(:,k,i),dveldx(:,k,i),rmatrix,ddenom)
          enddo
       !   print*,i,' exact derivs, dvxdy = ',dveldx(2,1,i)
       else
          print*,'WARNING: denominator collapsed in exact linear deriv = ',denom
          dveldx(:,:,i) = dveldx(:,:,i)*rho1i
       endif
       !read*
    endif

!
!--if evolving B instead of B/rho, add the extra term from the continuity eqn
!  (note that the dBevoldt term should be divided by rho)
!
    select case(imhd)
    case(11:19,21:) ! evolving B
       dBevoldt(:,i) = sqrtg(i)*dBevoldt(:,i) + Bevol(:,i)*rho1i*drhodt(i)
       if (idivBzero.ge.2) then
          !gradpsi(:,i) = gradpsi(:,i)*rho1i
          gradpsi(:,i) = gradpsi(:,i)*rhoi
          if (nsubsteps_divB.le.0) then
             dBevoldt(:,i) = dBevoldt(:,i) + gradpsi(:,i)
          endif
       endif
    case(10,20) ! remapped B/rho or remapped B
       dBevoldt(:,i) = 0.
    case(1:9) ! evolving B/rho
       if (iuse_exact_derivs.gt.0) then
          dBevoldt(:,i) = sqrtg(i)*dBevoldt(:,i)*rho1i
          !--add the B/rho dot grad v bit
          !if (any(dBevoldt(:,i).gt.0.)) then
          !print*,i,'dBevol/dt = ',dBevoldt(:,i)
          !endif
          do k=1,ndimV
             dBevoldt(k,i) = dBevoldt(k,i) + dot_product(Bevol(1:ndim,i),dveldx(1:ndim,k,i))
          enddo
          !if (any(dBevoldt(:,i).gt.0.)) then
          !print*,i,'dBevol/dt (exact) = ',dBevoldt(:,i)
          !read*
          !endif
       else
          dBevoldt(:,i) = sqrtg(i)*dBevoldt(:,i)*rho1i
          if (idivBzero.ge.2) then
             gradpsi(:,i) = gradpsi(:,i)*rho1i**2
          endif
       endif
    case(-1) ! vector potential evolution, crap gauge
       !
       !--add the v x B term
       !
       call cross_product3D(vel(:,i),Bfield(:,i),curlBi)
       dBevoldt(:,i) = dBevoldt(:,i)*rho1i + curlBi(:)
    case(-2) ! vector potential evolution, Axel gauge
       if (iuse_exact_derivs.gt.0) then
          do k=1,ndim
             dBevoldt(k,i) = -dot_product(Bevol(1:ndimV,i),dveldx(k,1:ndimV,i))
          enddo
       else
          dBevoldt(:,i) = dBevoldt(:,i)*rho1i
       endif
       !
       !--get v x Bext
       !
       call cross_product3D(vel(:,i),Bconst(:),curlBi)
       !--add v x Bext plus the existing term, which includes dissipation
       dBevoldt(:,i) = dBevoldt(:,i) + curlBi(:)
       !
       !--add dissipation for vector potential = -eta J
       !
       if (iav.ge.2) then
          curlBi(:) = curlB(:,i)
          !curlBi(:) = curlBsym(:,i)*rhoi
          !curlB(:,i) = curlBi(:)

          !curr2 = abs(dot_product(curlBsym(:,i)*rhoi,curlBsym(:,i)*rhoi))
          !curr2 = abs(dot_product(curlBi,curlBsym(:,i)*rhoi))
          curr2 = dot_product(curlBi,curlBi)
          etai = alpha(3,i)*sqrt(valfven2i)*hh(i)
          !etai = alpha(3,i)*vsig*hh(i)
          !etai = alpha(3,i)*2.*sqrt(curr2/rho(i))*hh(i)**2
          !etai = alpha(3,i)*sqrt(dot_product(graddivv(:,i),graddivv(:,i)))*hh(i)**2
          !etai = alpha(3,i)*(sqrt(valfven2i)*hh(i) + 2.*sqrt(curr2/rho(i))*hh(i)**2)
          !divB(i) = etai

          dBevoldt(:,i) = dBevoldt(:,i) - etai*curlBi(:)
          !if (curlBi(1).ne.0.) print*,i,curlB(1,i),curlBi(1)
          !print*,i,'diss = ',-etai,curlB(:,i)
          !
          !--add term to thermal energy equation. This is used to construct
          !  dendt for the entropy/total energy equations below
          !
          dudt(i) = dudt(i) + etai*curr2/rho(i)
       endif
    case(-3) ! Generalised Euler Potentials evolution
       dBevoldt(:,i) = 0.
    case default
       dBevoldt(:,i) = 0.
    end select
!
!--calculate resistive timestep (bootstrap onto force timestep)
!
    if (iresist.gt.0 .and. iresist.ne.2 .and. etamhd.gt.tiny(etamhd)) then
       if (iresist.eq.3) then
          etai = etafunc(x(1,i),etamhd)
       elseif (iresist > 0) then
          etai = etamhd
       endif
       dtforce = min(dtforce,hh(i)**2/etai)
    endif
!
!--if using the thermal energy equation, set the energy derivative
!  (note that dissipative terms are calculated in rates, but otherwise comes straight from cty)
!
    if (iener.eq.3) then
       ! could do this in principle but does not work with
       ! faniso modified by subtraction of Bconst
       if (damp.lt.tiny(0.)) dudt(i) = dudt(i) + pr(i)*rho1i**2*drhodt(i)
       dendt(i) = dot_product(vel(:,i),fprev(:,i)) + dudt(i) !&
         !       + 0.5*(dot_product(Bfield(:,i),Bfield(:,i))*rho1i**2) &
         !       + dot_product(Bfield(:,i),dBevoldt(:,i))*rho1i
    elseif (iener.eq.1) then ! entropy variable (just dissipative terms)
       if (damp.lt.tiny(0.)) dendt(i) = dendt(i) + (gamma-1.)/dens(i)**(gamma-1.)*dudt(i)
    elseif (iener.eq.4) then
       dendt(i) = (pr(i) + en(i))*rho1i*drhodt(i) + rhoi*dudt(i)
       if (iav.lt.0) stop 'iener=4 not compatible with Godunov SPH'
    elseif (iener.gt.0 .and. iav.ge.0 .and. iener.ne.10 .and. iener.ne.11 .and. &
           (idust.ne.1 .and. idust.ne.3 .and. idust.ne.4)) then
       if (damp.lt.tiny(0.)) dudt(i) = dudt(i) + pr(i)*rho1i**2*drhodt(i)
       dendt(i) = dudt(i)
    else
       dendt(i) = dudt(i)
    endif
    if (itype(i).eq.itypedust) dendt(i) = 0.
!
!--calculate time derivative of alpha (artificial dissipation coefficients)
!  see Morris and Monaghan (1997) and Price and Monaghan (2004c)
!
    daldt(:,i) = 0.
    if (any(iavlim.ne.0)) then
       tdecay1 = (avdecayconst*vsig)/hh(i)      ! 1/decay time (use vsig)
       !
       !--artificial viscosity parameter
       !
       select case(iavlim(1))
       case(1,2)
          source = max(drhodt(i)*rho1i,0.0)
          if (iavlim(1).eq.2) source = source*(2.0-alpha(1,i))
          daldt(1,i) = (alphamin - alpha(1,i))*tdecay1 + avfact*source
       case(3)
          graddivvmag = sqrt(dot_product(graddivv(:,i),graddivv(:,i)))
          !!print*,'graddivvmag = ',graddivvmag,max(drhodt(i)*rho1i,0.0)
          source = hh(i)*graddivvmag*(2.0-alpha(1,i))
          daldt(1,i) = (alphamin - alpha(1,i))*tdecay1 + avfact*source
       end select
       !
       !--artificial thermal conductivity parameter
       !
       if (iener.gt.0 .and. iavlim(2).gt.0) then
          !--this is using h*/sqrt(u)*(del^2 u) as the source
          if (uu(i).gt.epsilon(uu(i))) then
             sourceu = hh(i)*abs(del2u(i))/sqrt(uu(i))
             !sourceu = sqrt(abs(del2u(i))) !!
          else
             sourceu = 0.
          endif
          daldt(2,i) = (alphaumin - alpha(2,i))*tdecay1 + sourceu
       endif
       !
       !--artificial resistivity parameter
       !
       if (iavlim(3).ne.0 .and. imhd.ne.0) then
          !--calculate source term for the resistivity parameter
          sourceJ = SQRT(DOT_PRODUCT(curlB(:,i),curlB(:,i))*rho1i)
!          B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))
!          sourceJ = SQRT(vsig2*DOT_PRODUCT(curlB(:,i),curlB(:,i))/B2i)
          sourcedivB = 10.*abs(divB(i))*SQRT(rho1i)
          sourceB = MAX(sourceJ,sourcedivB)
          !sourceB = max(sqrt(drhodt(i)*rho1i*sourceJ),0.0)
          if (iavlim(3).eq.2) then
             sourceB = sourceB*(2.0-alpha(3,i))
          elseif (iavlim(3).eq.3) then
             source = max(drhodt(i)*rho1i,0.0)*(2.-alpha(3,i))
             !graddivv(:,i) = graddivv(:,i)/rho(i)
             !sourceJ = sqrt(dot_product(graddivv(:,i),graddivv(:,i)))
             !sourceB = sourceJ/(abs(drhodt(i)*rho1i) + sourceJ + tiny(0.))*sourceB*(2.-alpha(3,i))
             sourceB = sqrt(source*sourceB)
          endif
          daldt(3,i) = (alphaBmin - alpha(3,i))*tdecay1 + sourceB
       endif
    endif
!
!--calculate time derivative of divergence correction parameter psi
!
    select case(idivBzero)
       case(2:7)
          dpsidt(i) = -vsig2max*divB(i) - psidecayfact*psi(i)*vsigmax/hh(i)
          !if (abs(psi(i)).gt.0.) print*,' idivBzero',i,vsig2max,divB(i),psi(i)
       case DEFAULT
          dpsidt(i) = 0.
    end select

    if (idust.eq.1) then
       !
       !--DEBUGGING: CHECK ENERGY CONSERVATION
       !  BY ADDING TERMS (SHOULD GIVE ZERO)
       !
       esum = esum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - dot_product(vel(:,i),fexternal(:)) &
             + rhogasi*rhodusti(1)*rho1i**2*dot_product(deltav(:,i),ddeltavdt(:,i)) &
             + ((1. - 2.*dustfraci(1))*0.5*dot_product(deltav(:,i),deltav(:,i)) - uu(i))*ddustevoldt(1,i) &
             + rhogasi*rho1i*dudt(i))
    elseif (idust.eq.4) then
       esum = esum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - dot_product(vel(:,i),fexternal(:)) &
             - uu(i)*ddustevoldt(1,i) &
             + (1. - dustfrac(1,i))*dudt(i))
    endif
 enddo
 if ((idust.eq.1 .or. idust.eq.4) .and. abs(esum).gt.epsilon(esum) .and. (iener.ge.2)) then
    print*,' SUM (should be zero if conserving energy) = ',esum
    !read*
 endif
 if (idust.ne.0 .and. ratio > 1.) print "(a,g8.3,a)",' WARNING: max ts/dt = ',ratio,' approximation not valid'

 if (ivisc.gt.0) print*,' dEk/dt = ',esum

 if (sqrt(dot_product(fmean,fmean)).gt.1.e-8 .and. mod(nsteps,100).eq.0) print*,'WARNING: fmean = ',fmean(:)
!
!--calculate timestep constraint from the forces
!  dtforce is returned together with dtcourant to the main timestepping loop
!
 if (fhmax.lt.0.) then
    write(iprint,*) 'rates: fhmax <=0 :',fhmax
    call quit
 elseif (fhmax.gt.0.) then
    dtforce = min(dtforce,sqrt(1./fhmax))
 endif
 !!print*,'dtcourant = ',dtcourant,dtcourant2,0.2*dtcourant2
 !!dtcourant = 0.2*dtcourant2
!
!--set rates to zero on ghosts/fixed particles
!
 do i=1,ntotal      ! using ntotal just makes sure they are zero for ghosts
    if (itype(i).eq.itypebnd .or. itype(i).eq.itypebnddust .or. i.gt.npart) then
       force(:,i) = 0.0
       drhodt(i) = 0.0
       dhdt(i) = 0.0
       dudt(i) = 0.0
       dendt(i) = 0.0
       dBevoldt(:,i) = 0.0
       daldt(:,i) = 0.
       dpsidt(i) = 0.0
       fmag(:,i) = 0.0
       divB(i) = 0.0
       if (imhd.ge.0) curlB(:,i) = 0.0
       xsphterm(:,i) = 0.0
       gradpsi(:,i) = 0.0
    endif
 enddo

 if (allocated(listneigh)) deallocate(listneigh)
 if (allocated(phi)) deallocate(phi,del2u)
 if (trace) write(iprint,*) ' Exiting subroutine get_rates'
 if (itiming) then
    call cpu_time(t5)
    write(iprint,"(50('-'))")
    write(iprint,*) 'time for intro   = ',t2-t1,'s'
    write(iprint,*) 'time for main    = ',t3-t2,'s'
    write(iprint,*) 'time for gravity = ',t4-t3,'s'
    write(iprint,*) 'time for final   = ',t5-t4,'s'
    write(iprint,*) 'total rates time = ',t5-t1,'s'
    write(iprint,"(50('-'))")
 endif

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

contains
 real function slope_limiter(sl,sr,ilimiter) result(s)
  real, intent(in) :: sl,sr
  integer, intent(in) :: ilimiter

  s = 0.
  select case(ilimiter)
  case(6)
     if (sl*sr > 0.) s = 0.5*(sl + sr)
     ! van Albada
     !s = (sl + sl*sr)/(1. + sl*sr)
  case(5)
     ! UMIST
     if (sl*sr > 0.) s = sign(1.0,sl)*min(2.*abs(sr),&
                         (0.25*abs(sl) + 0.75*abs(sr)),&
                         (0.75*abs(sl) + 0.25*abs(sr)),2.*abs(sl))
  case(4)
     ! Superbee
     if (sl*sr > 0.) s = sign(1.0,sl)*max(min(abs(sr),2.*abs(sl)),min(2.*abs(sr),abs(sl)))
  case(3)
     ! minmod
     if (sl > 0. .and. sr > 0.) then
        s = min(abs(sl),abs(sr))
     elseif (sl < 0. .and. sr < 0.) then
        s = -min(abs(sl),abs(sr))
     endif
  case(2)
     ! Van Leer monotonised central (MC)
     if (sl*sr > 0.) s = sign(1.0,sl)*min(abs(0.5*(sl + sr)),2.*abs(sl),2.*abs(sr))
  case(1)
     ! van leer
     if (sl*sr > 0.) s = 2.*sl*sr/(sl + sr)
  case default
     ! no limiter
     s = 0.5*(sl + sr)
  end select

 end function slope_limiter

 subroutine reconstruct_dv(ndimV,ndim,vi,vj,dr,dx,dvdxi,dvdxj,projvstar,ilimiter,i,j)
  integer, intent(in)  :: ndimV,ndim,ilimiter
  real,    intent(in)  :: vi(ndimV),vj(ndimV),dr(ndimV),dx(ndim)
  real,    intent(in)  :: dvdxi(ndimV,ndimV),dvdxj(ndimV,ndimV)
  real,    intent(out) :: projvstar
  integer, intent(in), optional :: i,j
  real :: projvi,projvj,projv,si,sj,slope,projdvi(ndimV),projdvj(ndimV),absdx

  !if (ndim > 1) stop 'reconstruction not implemented for > 1D'
  projv = dot_product(vi,dr) - dot_product(vj,dr)
  !projvstar = projv - 0.5*dx(1)*dr(1)*(dvdxi(1,1) + dvdxj(1,1))

  do k=1,ndimV
     projdvi(k) = dot_product(dvdxi(k,:),dr)
     projdvj(k) = dot_product(dvdxj(k,:),dr)
  enddo
  si = dot_product(projdvi,dr)
  sj = dot_product(projdvj,dr)
  !si = dvdxj(1,1)
  !sj = dvdxj(1,1)
  absdx = sqrt(dot_product(dx,dx))
  !sep = dx(1)*dr(1)
  if (ilimiter > 0 .and. ilimiter < 10) then
     slope = slope_limiter(si,sj,ilimiter)
  else
     slope = 0.5*(si + sj)
  endif
  projvi = dot_product(vi,dr) - 0.5*absdx*slope
  projvj = dot_product(vj,dr) + 0.5*absdx*slope
  !projvstar = projvi - projvj

  projvstar = projv - absdx*slope
  !if (abs(absdx*slope) > abs(projv)) then
     !projvstar = sign(1.0,projv)*min(abs(projv),abs(projvstar))
     !print*,projvstar,projv,dot_product(vi,dr),projvi,dot_product(vj,dr),projvj
     !read*
  !endif
  !
  ! case 10 gives entropy limiter (projv and projvstar must have same sign)
  !
  if (ilimiter==10 .and. projvstar*projv < 0.) projvstar = projv !min(abs(projv),abs(projvstar))

 end subroutine reconstruct_dv

!----------------------------------------------
! Evaluate drag forces on a pair of particles
!----------------------------------------------
  subroutine drag_forces
    use kernels, only:interpolate_kerneldrag,interpolate_kernel
    use options, only:idrag_nature,Kdrag
    use dust,    only:get_tstop
    implicit none
    integer :: itypej
    real    :: dv2,vij,projv,projvstar,dragterm,dragterm_en
    real    :: spsoundgas,ts
    real    :: wabj,wab,hfacwabj,rho1j,gkerni,gkernj
    real, dimension(ndimV) :: drdrag
    real, parameter :: pi  = 4.*atan(1.)
!
!--differential velocities
!
    velj(:) = vel(:,j)
    dvel(:) = veli(:) - velj(:)
    dv2     = dot_product(dvel,dvel)
    if (rij.le.epsilon(rij)) then ! two particles at the same position
       if (dv2.le.epsilon(dv2)) then ! no differential velocity => no drag
          drdrag(1:ndim) = 0.
          return
       else ! Change dr so that the drag is colinear to the differential velocity
          vij       = sqrt(dv2)
          drdrag(:) = dvel(:)/vij
       endif
    else ! dr = drdrag
       drdrag(:) = dr(:)
    endif
!
!--calculate the kernel(s)
!
    if (idust.eq.2 .and. ibound(1).eq.5 .and. ndim <= 2) then
!--DIRTY HACK
       call interpolate_kernel(q2i,wabi,gkerni)
       call interpolate_kernel(q2j,wabj,gkernj)
    else
!--OK
       call interpolate_kerneldrag(q2i,wabi)
       call interpolate_kerneldrag(q2j,wabj)
    endif
    wabi     = wabi*hfacwabi
    hfacwabj = (1./hh(j)**ndim)
    wabj     = wabj*hfacwabj
!
!--get particle j properties
!
    itypej = itype(j)
    pmassj = pmass(j)
    rhoj   = rho(j)
    rho1j  = 1./rhoj
!
!--calculate the stopping time for use in the drag force
!
    projv = dot_product(dvel,drdrag)
    if (islope_limiter >= 0) then
       call reconstruct_dv(ndimV,ndim,vel(:,i),vel(:,j),dr,dx,gradB(:,:,i),gradB(:,:,j),projvstar,islope_limiter)
    else
       projvstar = projv
    endif
    if (itypei.eq.itypegas .or. itypei.eq.itypebnd) then
       spsoundgas = spsound(i)
       wab = wabi
       ts = get_tstop(idrag_nature,rhoi,rhoj,spsoundgas,Kdrag)
       h_on_csts_max = max(h_on_csts_max,hh(i)/(spsoundgas*ts))
    else
       if (itypej.ne.itypegas .and. itypej.ne.itypebnd) return
       spsoundgas = spsound(j)
       wab = wabj
       ts = get_tstop(idrag_nature,rhoj,rhoi,spsoundgas,Kdrag)
       h_on_csts_max = max(h_on_csts_max,hh(j)/(spsoundgas*ts))
    endif
    ts_min = min(ts_min,ts)
!
!--update the force and the energy
!
    if (idust.eq.2 .and. ibound(1).eq.5 .and. ndim <= 2) then
! HACK
       dragterm    = ndim*wab/((rhoi + rhoj)*ts)
       dragterm_en = dragterm*dot_product(dvel,dvel)
       forcei(:)   = forcei(:)  - dragterm*pmassj*dvel(:)
       force(:,j)  = force(:,j) + dragterm*pmassi*dvel(:)
    else
       dragterm    = ndim*wab/((rhoi + rhoj)*ts)*projvstar
       dragterm_en = dragterm*projv !dvdotr
       forcei(:)   = forcei(:)  - dragterm*pmassj*drdrag(:)
       force(:,j)  = force(:,j) + dragterm*pmassi*drdrag(:)
       if (dragterm_en < -1.) print*,' warning: dragterm_en = ',dragterm_en,i,j,projv,projvstar/projv
    endif
    if (itypei.eq.itypegas) then
       dudt(i)   = dudt(i) + pmassj*dragterm_en
    endif
    if (itypej.eq.itypegas) then
       dudt(j)   = dudt(j) + pmassi*dragterm_en
    endif

  end subroutine drag_forces

!--------------------------------------------------------------------------------------
! This is the interaction between the particle pairs
! Note that local variables are used from get_rates
!--------------------------------------------------------------------------------------
  subroutine rates_core
    use kernels, only:interpolate_kernel,interpolate_kernels
    use options, only:usenumdens
    use cons2prim, only:specialrelativity
    implicit none
    integer :: level
    real :: prstar,prstari,prstarj,vstar
    real :: dvsigdtc
    real :: hav,hav1,h21,q2
    real :: hfacwab,hfacwabj,hfacgrkern,hfacgrkernj
    real :: wabalti,wabaltj,wabalt
    real :: altrhoi,altrhoj,gammastar,vperp2
    real :: enthalpi,enthalpj,term,denom,term1
    real :: abs_dvcrossr,dvcrossr(ndimV)

   pmassj = pmass(j)
!
!--calculate the kernel(s)
!
   if (ikernav.eq.1) then
      hav = 0.5*(hi + hj)
      hav1 = 1./hav
      h21 = hav1*hav1
      hfacwab = hav1**ndim
      hfacgrkern = hfacwab*hav1
      q2 = rij2*h21
      call interpolate_kernel(q2,wab,grkern)
      wab = wab*hfacwab
      grkern = grkern*hfacgrkern
      grkerni = grkern
      grkernj = grkern
   else
      !  (using hi)
      call interpolate_kernels(q2i,wabi,grkerni,grkernalti,grgrkernalti)
      wabi = wabi*hfacwabi
      wabalti = wabi !wabalti*hfacwabi
      grkerni = grkerni*hfacgrkerni
      grkernalti = grkernalti*hfacgrkerni
      grgrkernalti = grgrkernalti*hfacgrkerni*hi1
      !  (using hj)
      hfacwabj = hj1**ndim
      hfacgrkernj = hfacwabj*hj1
      call interpolate_kernels(q2j,wabj,grkernj,grkernaltj,grgrkernaltj)
      wabj = wabj*hfacwabj
      wabaltj = wabj !wabaltj*hfacwabj
      grkernj = grkernj*hfacgrkernj
      grkernaltj = grkernaltj*hfacgrkernj
      grgrkernaltj = grgrkernaltj*hfacgrkernj*hj1
      !  (calculate average)
      wab = 0.5*(wabi + wabj)
      wabalt = 0.5*(wabalti + wabaltj)
      !  (grad h terms)
      if (ikernav.eq.3) then  ! if using grad h corrections
         if (usenumdens) then
         !--use these two lines for number density formulation
            grkerni = grkerni*(1. + gradhni*gradhi/pmassj)
            grkernj = grkernj*(1. + gradhn(j)*gradh(j)/pmassi)
         else
         !--use these two lines for usual formulation
            grkerni = grkerni*gradhi
            grkernj = grkernj*gradh(j)
         endif
         grkern = 0.5*(grkerni + grkernj)
      else  ! if not using grad h correction
         grkern = 0.5*(grkerni + grkernj)
         grkerni = grkern
         grkernj = grkern
      endif
   endif
!
!--define local copies of quantities
!
    velj(:) = vel(:,j)
    v2j = dot_product(velj(:),velj(:))
    dvel(:) = veli(:) - velj(:)
    dvdotr = dot_product(dvel,dr)
    call cross_product3D(dvel,dr,dvcrossr)
    abs_dvcrossr = sqrt(dot_product(dvcrossr,dvcrossr))
    projvi = dot_product(veli,dr)
    projvj = dot_product(velj,dr)
    rhoj = rho(j)
    rho1j = 1./rhoj
    rho2j = rhoj*rhoj
    rho21j = rho1j*rho1j
!!    rhoj5 = sqrt(rhoj)
    rhoij = rhoi*rhoj
    rhoav1 = 0.5*(rho1i + rho1j)   !2./(rhoi + rhoj)
    if (onef_dust) then
       dustfracj(:) = dustfrac(:,j)
       if (use_smoothed_rhodust) then
          rhodustj(:) = rhodust(:,j)
          rhogasj     = rhogas(j)
       else
          rhodustj(:)    = rhoj*dustfracj(:)
          rhogasj        = (1. - sum(dustfracj))*rhoj
       endif
       rhogrhodonrhoj(:) = rhogasj*rhodustj(:)*rho1j
       deltavj(:)     = deltav(:,j)
       deltav2j       = dot_product(deltavj,deltavj)
       vgasi(:)       = veli(:) - sum(dustfraci)*deltavi(:)
       vgasj(:)       = velj(:) - sum(dustfracj)*deltavj(:)
       dvgas(:)       = vgasi(:) - vgasj(:)
       projdvgas      = dot_product(dvgas,dr)
       projdeltavi    = dot_product(deltavi,dr)
       projdeltavj    = dot_product(deltavj,dr)
    else
       projdvgas      = dvdotr
       deltav2j       = 0.
       rhogasj        = rhoj
    endif
    prj = max(pr(j) - pext,0.)
    prnetj = prj - pequil(iexternal_force,x(:,j),rhoj)
    Prho2j = prj*rho21j
    spsoundj = spsound(j)
    uuj = uu(j)

    phii_on_phij = phii/phi(j)
    phij_on_phii = phi(j)*phii1
    sqrtgj = sqrtg(j)
    !-- mhd definitions --
    if (imhd.ne.0) then
       Bj(:) = Bfield(:,j)
       Brhoj(:) = Bj(:)*rho1j
       dB(:) = Bi(:) - Bj(:)
       projBi = dot_product(Bi,dr)
       projBj = dot_product(Bj,dr)
       projdB = dot_product(dB,dr)
       projBrhoi = dot_product(Brhoi,dr)
       projBrhoj = dot_product(Brhoj,dr)

       if (trim(geom).ne.'cartes') then
          call metric_diag(x(:,j),gdiagj(:),sqrtgj,ndim,ndimV,geom)
          B2j = dot_product_gr(Bj,Bj,gdiagj)
          valfven2j = B2j/dens(j)
       else
          B2j = dot_product(Bj,Bj)
          valfven2j = B2j*rho1j
       endif
       Brho2j = B2j*rho21j
       projBconst = dot_product(Bconst,dr)
       if (imhd.lt.0) then
          dBevol(:) = Bevoli(:) - Bevol(:,j)
       endif
       if (iresist.eq.3) then
          etaj = etafunc(x(1,j),etamhd)
       elseif (iresist > 0) then
          etaj = etamhd
       endif
    endif
    forcej(:) = 0.
    fextraj(:) = 0.

    !--maximum velocity for timestep control
    !            vmag = SQRT(DOT_PRODUCT(dvel,dvel))
    !            vsigdtc = vmag + spsoundi + spsoundj        &
    !                         + sqrt(valfven2i) + sqrt(valfven2j)
    vsig = 0.
    vsigav = 0.
!----------------------------------------------------------------------------
!  calculate signal velocity (this is used for timestep control
!  and also in the artificial viscosity)
!----------------------------------------------------------------------------
    if (specialrelativity) then
       !
       ! special relativistic signal velocities
       !
       enthalpi = 1. + uui + pri/dens(i)
       enthalpj = 1. + uuj + prj/dens(j)
       spsoundi = sqrt(gamma*pri/(dens(i)*enthalpi))
       spsoundj = sqrt(gamma*prj/(dens(j)*enthalpj))
       select case(iav)
       case(5)
          vsigi = max((veli(1) + spsoundi)/(1. + veli(1)*spsoundi), &
                      (veli(1) - spsoundi)/(1. - veli(1)*spsoundi))
          vsigj = max((velj(1) + spsoundj)/(1. + velj(1)*spsoundj), &
                      (velj(1) - spsoundj)/(1. - velj(1)*spsoundj))
          vsig = max(vsigi,vsigj,0.)
       case(4)
          !
          ! just using the speed of light...
          !
          vsig = 1.
       case(3)
          !
          ! this is vsig(1) from Chow & Monaghan '97 (equation 5.7)
          !
          vstar= abs((projvi-projvj)/(1.0-(projvi*projvj))) ! (equation 4.20)
          vsigi=(spsoundi+vstar)/(1.0+(spsoundi*vstar))
          vsigj=(spsoundj+vstar)/(1.0+(spsoundj*vstar))
          !pmomstar
          vsig=min(vsigi+vsigj+vstar,1.)

       case(2)
          !
          ! vsig(2) from Chow & Monaghan '97 (equation 5.9)
          !
          vstar= abs((projvi-projvj)/(1.0-(projvi*projvj))) ! (equation 4.20)
          gammastar=1.0/(sqrt(1.0-vstar**2))

          !vsig=1
          vsigi=sqrt(((spsoundi**2)+beta*(gammastar-1.0))/(1+((spsoundi**2)+((beta*(gammastar-1))/(gamma-1)))))
          vsigj=sqrt(((spsoundj**2)+beta*(gammastar-1.0))/(1+((spsoundj**2)+((beta*(gammastar-1))/(gamma-1)))))
          vsig=abs(((vsigi+projvi)/(1.+vsigi*projvi))-((vsigj-projvj)/(1.-vsigj*projvj)))

       case default

          !
          ! vsig from the maximum eigenvalue, (e.g. Font et al., Rosswog 2010, Joe 2001)
          !

          vperp2 = v2i - projvi**2
          denom  = (1. - v2i*spsoundi**2)
          term   = spsoundi*sqrt((1.-v2i)*(1.-projvi**2 - spsoundi**2*vperp2))
          term1  = projvi*(1.-spsoundi**2)
          vsigi  = max(term1 + term,abs(term1 - term))/denom

          vperp2 = v2j - projvj**2
          denom  = (1. - v2j*spsoundj**2)
          term   = spsoundj*sqrt((1.-v2j)*(1.-projvj**2 - spsoundj**2*vperp2))
          term1  = projvj*(1. - spsoundj**2)
          vsigj  = max(term1 + term,abs(term1 - term))/denom

          vsig = max(vsigi,vsigj)
          !print*,' vsig = ',vsigi,vsigj,spsoundi,spsoundj

       end select
       if (vsig.gt.1) stop 'error: vsig > 1'
       vsigu = vsig
       vsigdtc = vsig

    else
       !
       !--max signal velocity (Joe's)
       !
   !    vsigi = 0.5*(sqrt(spsoundi**2 + valfven2i - 2.*spsoundi*projBi/rhoi5) &
   !                +sqrt(spsoundi**2 + valfven2i + 2.*spsoundi*projBi/rhoi5))
   !    vsigj = 0.5*(sqrt(spsoundj**2 + valfven2j - 2.*spsoundj*projBj/rhoj5) &
   !                +sqrt(spsoundj**2 + valfven2j + 2.*spsoundj*projBj/rhoj5))

       !
       !--max signal velocity (my version)
       !
       if (imhd.ne.0) then
          vsig2i = spsoundi**2 + valfven2i
          vsig2j = spsoundj**2 + valfven2j
          vsigproji = vsig2i**2 - 4.*(spsoundi*(projBi))**2*rho1i
          vsigprojj = vsig2j**2 - 4.*(spsoundj*(projBj))**2*rho1j
          if (vsigproji.lt.0.) then
             write(iprint,*) ' rates: i=',i,' vsig det < 0 ',vsigproji,vsig2i**2,4*(spsoundi*projBi)**2*rho1i
             call quit
          elseif (vsigprojj.lt.0.) then
             write(iprint,*) ' rates: j=',j,' vsig det < 0 ',vsigprojj,vsig2j**2,4*(spsoundj*projBj)**2*rho1j
             print*,j,'rho,dens = ',rho(j),dens(j)
             call quit
          endif
          vsigi = sqrt(0.5*(vsig2i + sqrt(vsigproji)))
          vsigj = sqrt(0.5*(vsig2j + sqrt(vsigprojj)))
          if (iavlim(3).ne.2) then
             vsigB = norm2(dvel) ! abs_dvcrossr + abs(dvdotr)
          else
             vsigB = 0.5*(vsigi + vsigj) + abs(dvdotr)
          endif

          !vsigB = sqrt(dot_product(dvel - dvdotr,dvel - dvdotr))
          !vsigB = 0.5*(sqrt(valfven2i) + sqrt(valfven2j))
       elseif (iquantum > 0) then
          vsigi = sqrt(spsoundi**2 + 1.0 / (4.0 * hi**2))
          vsigj = sqrt(spsoundj**2 + 1.0 / (4.0 * hj**2))
          vsigB = 0.
          !vsigi = spsoundi
          !vsigj = spsoundj
       else
          vsigi = spsoundi
          vsigj = spsoundj
          vsigB = 0.
       endif

       vsig = 0.5*(max(vsigi + vsigj - beta*dvdotr,0.0)) ! also used where dvdotr>0 in MHD
       !vsig = sqrt(vsigi**2 + beta*dvdotr**2) + sqrt(vsigj**2 + beta*dvdotr**2) - dvdotr
       !vsig = max(4.*vsigi*vsigj/(vsigi + vsigj) - beta*dvdotr,0.)
       !vsig = 0.5*(vsigi + vsigj + beta*abs(dvdotr))
       !vsigB = 0.5*max(vsigi + vsigj - 4.0*dvdotr,0.0) !!!*(1./(1.+exp(1000.*dvdotr/vsig)))
       !vsigB = max(-dvdotr,0.0) !!!*(1./(1.+exp(1000.*dvdotr/vsig)))

       vsigu = sqrt(abs(prneti-prnetj)*rhoav1)
       !vsigu = sqrt(0.5*(psi(i) + psi(j)))
       !vsigu = abs(dvdotr)
       !vsigu = sqrt(abs(uui - uuj))

       ! vsigdtc is the signal velocity used in the timestep control
       vsigdtc = max(vsig,0.5*(vsigi + vsigj + beta*abs(dvdotr)),vsigB)
       if (idust.eq.1) then
          vsigdtc = vsigdtc + sqrt(deltav2i + deltav2j)
       endif

    endif

    if (itype(i).eq.itypedust) then
       vsig  = 0.
       vsigu = 0.
    else
       dvsigdtc = 1./vsigdtc
       vsigmax = max(vsigmax,vsigdtc)
    !
    !--time step control (courant and viscous)
    !
       if (vsigdtc.gt.zero) dtcourant = min(dtcourant,min(hi*dvsigdtc,hj*dvsigdtc))
    endif
!----------------------------------------------------------------------------
!  artificial dissipation terms
!----------------------------------------------------------------------------

    if (iav.gt.0) then
       if (specialrelativity) then
          call artificial_dissipation_sr
       elseif (idust.eq.1) then
          call artificial_dissipation_dust
       elseif (idust.eq.3 .or. idust.eq.4) then
          call artificial_dissipation_dust_diffusion
       elseif (iav.eq.3) then
          call artificial_dissipation_phantom
       else
          call artificial_dissipation
       endif
    endif
    if (vsigav.gt.zero) dtav = min(dtav,min(hi/vsigav,hj/vsigav))

!------------------------------------------------------------------------
!  Physical viscosity terms
!------------------------------------------------------------------------
    if (ivisc.gt.0) call physical_viscosity
    !if (vsigav.gt.zero) dtav = min(dtav,min(hi/vsigav,hj/vsigav))

!----------------------------------------------------------------------------
!  pressure term (generalised form)
!----------------------------------------------------------------------------
    if (iprterm.ge.0) then
       if (iav.lt.0) then
          !!if (abs(pri-prj).gt.1.e-2) then
          if (iav.eq.-2) then
             !--lanzafame (2009)-style
             if (dvdotr.lt.0.) then
                term = atan(-5.*dvdotr/spsoundi)/1.57077
                prstari = pri*(1.-alphai*term*dvdotr/spsoundi)**2
                term = atan(-5.*dvdotr/spsoundj)/1.57077
                prstarj = prj*(1.-alpha(1,j)*term*dvdotr/spsoundj)**2
             else
                prstari = pri
                prstarj = prj
             endif
          else
             call riemannsolver(gamma,pri,prj,-projvi,-projvj, &
               rhoi,rhoj,prstar,vstar)
             prstari = prstar
             prstarj = prstar
          endif
!             if (prstar.gt.1.00001*max(pri,prj) .or. prstar.lt.0.99999*min(pri,prj)) then
!                print*,'error, pstar = ',prstar,'i,j = ',pri,prj,projvi,projvj
!             endif
             !prstar = 0.5*(pri + prj)
          prterm = prstari*phii_on_phij*rho21i*sqrtgi*grkerni &
                 + prstarj*phij_on_phii*rho21j*sqrtgj*grkernj
       else
          prterm = phii_on_phij*Prho2i*sqrtgi*grkerni &
                 + phij_on_phii*Prho2j*sqrtgj*grkernj
          if (ibiascorrection.ne.0) then
             prterm = prterm + ((1. - sqrtgi)*pri/rhoalt(i)**2*grkernalti &
                              + (1. - sqrtgj)*prj/rhoalt(j)**2*grkernaltj)
          endif

!          print *, 'pr ', phii_on_phij,Prho2i,sqrtgi,grkerni, &
!                 phij_on_phii, Prho2j, sqrtgj,grkernj
          if (iprterm.eq.10) then
!
!--Ritchie and Thomas force eqn.
!
             !altrhoi = pri/((gamma-1.)*uui)
             !altrhoj = prj/((gamma-1.)*uuj)
             prterm = uui*uuj*(gamma-1.)**2*(1./pri*grkerni &
                                           + 1./prj*grkernj)
          elseif (iprterm.eq.12) then
!
!--volume correction from uhat
!
             prterm = psi(i)*Prho2i*sqrtgi*grkerni &
                    + psi(j)*Prho2j*sqrtgj*grkernj
          endif
       endif
       !
       !--add pressure terms to force
       !
!       print *, 'pres ', prterm, dr(1), dr(2)
       forcei(:) = forcei(:) - pmassj*prterm*dr(:) !- pmassj*wab*rij
       forcej(:) = forcej(:) + pmassi*prterm*dr(:) !+ pmassi*wab*rij
    elseif (iprterm.eq.-2) then
       !--use Pa-Pb pressure force
       forcei(:) = forcei(:) + pmassj*(pri - prj)*rho21j*dr(:)*grkern
       forcej(:) = forcej(:) + pmassi*(pri - prj)*rho21i*dr(:)*grkern
    endif
!------------------------------------------------------------------------
!  Lorentz force and time derivative of B terms
!------------------------------------------------------------------------

    if (imhd.ne.0) call mhd_terms

!------------------------------------------------------------------------
!  QSPH terms
!------------------------------------------------------------------------

    if (iquantum.ne.0) call quantum_terms

!------------------------------------------------------------------------
!  Time derivative of dust to gas ratio and deltav for one fluid dust
!------------------------------------------------------------------------

    if (idust.eq.1) then
       call dust_derivs
    elseif (idust.eq.4) then
       call dust_derivs_diffusion
    endif
!
!   Add contributions to j from mhd terms and dissipation here
!   to avoid repeated memory access
!
!   here forcej is the gas-only forces, fextraj are forces that act on the total fluid
    force(:,j) = force(:,j) + fextraj(:) + forcej(:)

    if (ind_timesteps.ne.0) then
       level = min(ibin(i),ibin(j))
       !print*,' level = ',level
       force_bins(:,j,level+1) = force_bins(:,j,level+1) + fextraj(:) + forcej(:)
       force_bins(:,i,level+1) = force_bins(:,i,level+1) - pmassj/pmassi*(fextraj(:) + forcej(:))
    endif

!------------------------------------------------------------------------
!  total energy equation (thermal energy equation terms calculated
!                         outside loop and in artificial_dissipation)
!------------------------------------------------------------------------

    if (iener.eq.3 .and. .false.) then
       call energy_equation
    elseif (iener.gt.0 .and. iav.lt.0) then
       dudt(i) = dudt(i) + pmassj*prstari*(rho21i)*dvdotr*grkerni
       dudt(j) = dudt(j) + pmassi*prstarj*(rho21j)*dvdotr*grkernj
    elseif (iener.eq.10) then
       !
       !--Ritchie & Thomas energy equation (their eqn. 25)
       !
       altrhoi = pri/((gamma-1.)*uui)
       altrhoj = prj/((gamma-1.)*uuj)
       dudt(i) = dudt(i) + pmassj*(gamma-1.)*uuj/altrhoi*dvdotr*grkerni
       dudt(j) = dudt(j) + pmassi*(gamma-1.)*uui/altrhoj*dvdotr*grkernj
    elseif (iener.eq.11) then
       !
       !--usual form, not directly related to drho/dt
       !
       dudt(i) = dudt(i) + Prho2i*pmassj*dvdotr*grkerni
       dudt(j) = dudt(j) + Prho2j*pmassi*dvdotr*grkernj
    endif

!------------------------------------------------------------------------
!  grad u term for dissipation switch
!------------------------------------------------------------------------

    if (iav.gt.0) then
       if(iavlim(2).gt.0) then
          graduterm = (uu(i)-uu(j))/rij
          del2u(i) = del2u(i) + pmassj*rho1j*graduterm*grkerni
          del2u(j) = del2u(j) - pmassi*rho1i*graduterm*grkernj
       endif
       if (iavlim(1).eq.3) then
          !!graddivvterm =
          graddivv(:,i) = graddivv(:,i) + pmassj*rho1j/rij*dvdotr*grkerni*dr(:)
          graddivv(:,j) = graddivv(:,j) - pmassi*rho1i/rij*dvdotr*grkernj*dr(:)
       else
          !
          !--calculate curl v
          !
          graddivv(:,i) = graddivv(:,i) + pmassj*(dvel(:) - dvdotr)*grkerni
          graddivv(:,j) = graddivv(:,j) + pmassi*(-dvel(:) - dvdotr)*grkernj
       endif
    endif

!----------------------------------------------------------------------------
!  XSPH term for moving the particles
!----------------------------------------------------------------------------
    if (ixsph.eq.1) then
       xsphterm(:,i) = xsphterm(:,i) - pmassj*dvel(:)*rhoav1*wabalt
       xsphterm(:,j) = xsphterm(:,j) + pmassi*dvel(:)*rhoav1*wabalt
!       dudt(i) = dudt(i) - 0.2*pmassj*(uu(i)-uu(j))*rhoav1*grkernalti
!       dudt(j) = dudt(j) + 0.2*pmassi*(uu(i)-uu(j))*rhoav1*grkernaltj
    endif
!
!  Continuity equation if density not done by summation
!
    if (icty.ge.1) then
       drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni
       drhodt(j) = drhodt(j) + pmassi*dvdotr*grkernj
    endif

!----------------------------------------------------------------------------
!  grad v
!----------------------------------------------------------------------------
    if (allocated(dveldx)) then
       do k=1,ndimV
          dveldx(:,k,i) = dveldx(:,k,i) - pmassj*dvel(k)*dr(:)*grkerni
          dveldx(:,k,j) = dveldx(:,k,j) - pmassi*dvel(k)*dr(:)*grkernj
       enddo
       do k=1,ndxdx
          dxdx(k,i) = dxdx(k,i) - pmass(j)*(dx(idxdx(k)))*dr(jdxdx(k))*grkerni
          dxdx(k,j) = dxdx(k,j) - pmass(i)*(dx(idxdx(k)))*dr(jdxdx(k))*grkernj
       enddo
    endif

    return
  end subroutine rates_core

!--------------------------------------------------------------------------------------
! These are the artificial viscosity, thermal conduction and resistivity terms
! Change this to change the artificial viscosity algorithm
! Inherits the values of local variables from rates
!
! This version corresponds to the MHD dissipative terms
! described in Price & Monaghan (2004a), MNRAS 348, 137
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation
    implicit none
    real, dimension(ndimB) :: Bvisc,dBdtvisc
    real :: visc,alphaav,alphaB,alphau
    real :: v2i,v2j,B2i,B2j
    real :: qdiff
    real :: vissv,vissB,vissu
    real :: term,dpmomdotr
    real :: termB,termu,termv
    !
    !--definitions
    !
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    alphaB = 0.5*(alphaBi + alpha(3,j))
    vsigav = max(alphaav,alphau,alphaB)*vsig
    !!rhoav1 = 2./(rhoi + rhoj)
    if (geom(1:4).ne.'cart') then
       dpmomdotr = abs(dot_product(pmom(:,i)-pmom(:,j),dr(:)))
    else
       dpmomdotr = -dvdotr
    endif
    !--used for viscosity
    term = vsig*rhoav1*grkern
    termv = term

    !--used for thermal conductivity
    termu = vsigu*rhoav1*grkern
    !termu = vsig*rho1i*rho1j*grkern
    !termu = vsig*rhoav1*grkern

    !--used for resistivity
    termB = vsigB*rhoav1*grkern
!    termB = vsigi*alphaBi*vsigj*alpha(3,j)/ &
!                (vsigi*alphaBi + vsigj*alpha(3,j))*rhoav1*grkern
!    termB = vsigi*alphaBi*vsigj*alpha(3,j)*grkerni*grkernj/ &
!                (rhoi*vsigi*alphaBi*grkerni + rhoj*vsigj*alpha(3,j)*grkernj)
    !termB = vsigB*rhoav1*grkern


    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------

    if (dvdotr.lt.0 .and. iav.le.3) then
       visc = alphaav*termv*dpmomdotr     ! viss=abs(dvdotr) defined in rates
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
    elseif (iav.eq.4) then ! using total energy, for approaching and receding
       visc = alphaav*termv
       !print*,'visc = ',i,j,vsig*alphaav*hh(i)
       !visc = visc*0.5*(hh(i) + hh(j))/rij  ! use this line to multiply viscosity by h/r
       forcei(:) = forcei(:) + pmassj*visc*dvel(:)
       forcej(:) = forcej(:) - pmassi*visc*dvel(:)
    endif

    !--------------------------------------------
    !  resistivity in induction equation
    !  (applied everywhere)
    !--------------------------------------------
    if (imhd.ne.0) then
       if (imhd.gt.0 .and. imhd.ne.10) then
          if (iav.ge.2) then
             Bvisc(:) = dB(:)*rhoav1
          else
             Bvisc(:) = (dB(:) - dr(:)*projdB)*rhoav1
          endif
          dBdtvisc(:) = alphaB*termB*Bvisc(:)

          !
          !--add to d(B/rho)/dt (converted to dB/dt later if required)
          !
          dBevoldti(:) = dBevoldti(:) + rhoi*pmassj*dBdtvisc(:)
          dBevoldt(:,j) = dBevoldt(:,j) - rhoj*pmassi*dBdtvisc(:)

       elseif (iav.eq.1) then !--wrong vector potential resistivity
          dBdtvisc(:) = alphaB*termB*dBevol(:)
          !
          !--add to dA/dt (note: do not multiply by rho here)
          !
          dBevoldti(:) = dBevoldti(:) + pmassj*dBdtvisc(:)
          dBevoldt(:,j) = dBevoldt(:,j) - pmassi*dBdtvisc(:)
       endif
    endif

    !--------------------------------------------------
    !  dissipation terms in energy equation
    !  (viscosity + resistivity + thermal conduction)
    !
    !   total energy equation
    !--------------------------------------------------
    if (iener.eq.3) then

       qdiff = 0.
       !
       !  kinetic energy terms - applied only when particles approaching
       !
       if (dvdotr.lt.0 .and. iav.le.3) then
          v2i = dot_product(veli,dr)**2      ! energy along line
          v2j = dot_product(velj,dr)**2      ! of sight
          qdiff = qdiff + term*alphaav*0.5*(v2i-v2j)
       elseif (iav.eq.4) then
          v2i = dot_product(veli,veli)      ! total energy
          v2j = dot_product(velj,velj)
          qdiff = qdiff + term*alphaav*0.5*(v2i-v2j)
       endif
       !
       !  thermal energy terms - applied everywhere
       !
       qdiff = qdiff + alphau*termu*(uu(i)-uu(j))
       !
       !  magnetic energy terms - applied everywhere
       !
       if (imhd.gt.0 .and. imhd.ne.10) then
          if (iav.ge.2) then
             B2i = dot_product(Bi,Bi) ! total magnetic energy
             B2j = dot_product(Bj,Bj)
          else
             B2i = (dot_product(Bi,Bi) - dot_product(Bi,dr)**2) ! magnetic energy
             B2j = (dot_product(Bj,Bj) - dot_product(Bj,dr)**2) ! along line of sight
          endif
          qdiff = qdiff + alphaB*termB*0.5*(B2i-B2j)*rhoav1
       elseif (imhd.lt.0) then
          stop 'mhd dissipation not implemented with total energy equation for vector potential'
       endif
       !
       !  add to total energy equation
       !
       dendt(i) = dendt(i) + pmassj*qdiff
       dendt(j) = dendt(j) - pmassi*qdiff

    !--------------------------------------------------
    !   thermal energy equation
    !--------------------------------------------------
    elseif (iener.gt.0) then
       !
       !  kinetic energy terms
       !
       if (dvdotr.lt.0 .and. iav.le.3) then
          vissv = -alphaav*0.5*(dot_product(veli,dr) - dot_product(velj,dr))**2
       elseif (iav.ge.4) then
          vissv = -alphaav*0.5*dot_product(dvel,dvel)
       else
          vissv = 0.
       endif
       !
       !  thermal energy terms
       !
       if (iener.eq.1) then
          vissu = 0.
       else
          !vissu = alphau*(pri - prj)
          vissu = alphau*(uu(i) - uu(j))
       endif
       !
       !  add magnetic energy term - applied everywhere
       !
       if (imhd.gt.0 .and. imhd.ne.10) then
          if (iav.ge.2) then
             vissB = -alphaB*0.5*(dot_product(dB,dB))*rhoav1
          else
             vissB = -alphaB*0.5*(dot_product(dB,dB)-projdB**2)*rhoav1
          endif
       elseif (imhd.lt.0 .and. iav.eq.1 .and. imhd.ne.-3) then
          !vissB = -alphaB*0.5*(dot_product(curlB(:,i)-curlB(:,j),dBevol(:)))*rhoav1
          vissB = -alphaB*0.5*(dot_product(dB,dB))*rhoav1
       else
          vissB = 0.
       endif
       !
       !  add to thermal energy equation
       !
       if (damp.lt.tiny(0.)) then
          dudt(i) = dudt(i) + pmassj*(term*(vissv) + termu*vissu + termB*(vissB))
          dudt(j) = dudt(j) + pmassi*(term*(vissv) - termu*vissu + termB*(vissB))

          !--entropy dissipation
          if (iener.eq.1) then
             if (alphau.gt.0.) stop 'thermal conduction in entropy equation is wrong'
             vissu = alphau*(en(i)-en(j))
             dendt(i) = dendt(i) + pmassj*(termu*(vissu))
             dendt(j) = dendt(j) + pmassi*(termu*(-vissu))
          endif
       endif

       !if (icty.eq.1) then
       !   vissrho = alphaB*vsig*grkern*(rhoi - rhoj)*rhoav1 !!/sqrt(rhoi*rhoj)
       !   drhodt(i) = drhodt(i) + pmassj*(vissrho)
       !   drhodt(j) = drhodt(j) + pmassi*(-vissrho)
       !endif
    endif

    return
  end subroutine artificial_dissipation

!--------------------------------------------------------------------------------------
! These are the artificial viscosity, thermal conduction and resistivity terms
! Change this to change the artificial viscosity algorithm
! Inherits the values of local variables from rates
!
! This version corresponds to av terms as implemented in phantom
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation_phantom
    implicit none
    real :: vsi,vsj,qi,qj,visc,du,cfaci,cfacj,diffu
    real :: dudti,dudtj,projv,prstar,vstar

    dudti = 0.
    dudtj = 0.

    if (dvdotr < 0. .or. islope_limiter >= 0) then
       !
       ! artificial viscosity, as in Phantom
       !
       vsi = max(alphai*spsoundi - beta*dvdotr,0.)
       vsj = max(alpha(1,j)*spsoundj - beta*dvdotr,0.)

       projv = dvdotr
       if (islope_limiter >= 0 .and. .false.) then
          call reconstruct_dv(ndimV,ndim,vel(:,i),vel(:,j),dr,dx,gradB(:,:,i),gradB(:,:,j),projv,islope_limiter,i,j)
!          call reconstruct_dv(ndimV,ndim,vel(:,i),vel(:,j),dr,dx,-drhodt(i)*rho1i,-drhodt(j)*rho1j,projv,islope_limiter,i,j)
          !if (-dvdotr > 0.1) print*,dvdotr,projv
       endif
       !call riemannsolver(gamma,prneti,prnetj,-projvi,-projvj, &
       !                   rhoi,rhoj,prstar,vstar)
       !qi = prstar - prneti
       !qj = prstar - prnetj
       qi = -0.5*rhoi*vsi*projv
       qj = -0.5*rhoj*vsj*projv

       visc = (qi*rho21i*grkerni + qj*rho21j*grkernj)
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
       !
       !  add to thermal energy equation
       !
       if (damp.lt.tiny(0.)) then
          dudti = qi*rho21i*pmassj*dvdotr*grkerni
          dudtj = qj*rho21j*pmassi*dvdotr*grkernj
       endif
    endif

    if (damp < tiny(damp)) then
       !
       ! artificial thermal conductivity, as in PL15
       !
       !call riemannsolver(gamma,prneti,prnetj,projvi,projvj, &
       !                   rhoi,rhoj,prstar,vstar)
       !call hllc_solver(gamma,prneti,prnetj,projvi,projvj,rhoi,rhoj,prstar,vstar)
       !vsigu = max(vstar,0.) !vsigu = max(vstar,0.) !abs(vstar)
       du = uu(i) - uu(j)
       cfaci = 0.5*alphaui*rhoi*vsigu*du
       cfacj = 0.5*alpha(2,j)*rhoj*vsigu*du
       diffu = cfaci*grkerni*rho1i**2 + cfacj*grkernj*rho1j**2
       dudt(i) = dudt(i) + dudti + pmassj*diffu
       dudt(j) = dudt(j) + dudtj - pmassi*diffu
    endif

  end subroutine artificial_dissipation_phantom

!--------------------------------------------------------------------------------------
! Artificial dissipation for one-fluid dust:
!  iav = 1 => gas-only dissipation, as in submitted LP14 paper
!  iav = 2 => M97-inspired dissipation, with heating from deltav term
!  iav = 3 => dissipation in du/dt is 0.5*(vab.rab)^2
!  iav = 4 => gas-only dissipation but without projecting velocities along line of sight
!
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation_dust
    implicit none
    real :: visc,alphaav,alphau,alphaB
    real :: vissv,vissu,vissdv,vsigdv
    real :: term,dpmomdotr,faci,facj
    real :: termu,termv,termdv,vsigeps,diffeps

    real :: dustfracav,projdvgasav,projddeltav,projepsdeltav
    real, dimension(ndimV) :: ddeltav
    logical :: allpairs
    !
    !--definitions
    !
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    alphaB = 0.5*(alphaBi + alpha(3,j))
    vsigav = max(alphaav,alphau)*vsig

    dustfracav = 0.5*(dustfraci(1) + dustfracj(1))
    projdvgasav = dvdotr - dustfracav*(projdeltavi - projdeltavj)
    ddeltav     = deltavi(:) - deltavj(:)
    projddeltav = dot_product(ddeltav,dr)
    projepsdeltav = dustfraci(1)*projdeltavi - dustfracj(1)*projdeltavj

    if (iav.eq.3) then
       dpmomdotr = projdvgasav
    elseif (iav.eq.2) then
       dpmomdotr = projdvgas
    else
       dpmomdotr = dvdotr
    endif
    !--used for viscosity
    term = vsig*rhoav1*grkern
    termv = term

    allpairs = .false.
    !allpairs = .true.
    if (iav.eq.4 .or. iav.eq.1) allpairs = .true.

    if (iav.eq.2 .or. iav.eq.3) then
       termv = termv*(1. - dustfracav)
    elseif (iav.eq.4) then
       termv = termv !*dustfracav*(1. - dustfracav)
       if (iav.gt.4) stop 'av terms not implemented for iav>4 and one-fluid dust'
    endif

    !--used for thermal conductivity
    termu = vsigu*rhoav1*grkern*(1. - dustfracav)  ! multiply by (1-dustfrac)

    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------
    vsigdv = vsig !0.5*sqrt(projdeltavi**2 + projdeltavj**2)

    if (projdvgas.lt.0 .or. allpairs) then
       visc = alphaav*termv*dpmomdotr     ! viss=abs(dvdotr) defined in rates
       if (iav.eq.4) then
          ! just to gas
          termdv = vsigdv*rhoav1*grkern
          fextrai(:) = fextrai(:) + pmassj*alphaav*(termv*dvdotr - termdv*projepsdeltav)*dr(:)
          fextraj(:) = fextraj(:) - pmassi*alphaav*(termv*dvdotr - termdv*projepsdeltav)*dr(:)
       elseif (iav.eq.1 .or. iav.eq.3) then
          ! here viscosity applies to the whole fluid, not just gas component
          fextrai(:) = fextrai(:) + pmassj*visc*dr(:)
          fextraj(:) = fextraj(:) - pmassi*visc*dr(:)
       else
          forcei(:) = forcei(:) + pmassj*visc*dr(:)
          forcej(:) = forcej(:) - pmassi*visc*dr(:)
       endif
    endif

    !--------------------------------------------------
    !  dissipation terms in energy equation
    !--------------------------------------------------
    if (iener.eq.3) then
       stop 'onefluid dust dissipation not implemented in total energy equation'
    !--------------------------------------------------
    !   thermal energy equation
    !--------------------------------------------------
    elseif (iener.gt.0) then
       !
       !  kinetic energy terms
       !
       if (projdvgas.lt.0 .or. allpairs) then
          if (iav.eq.4) then
             vissv = -alphaav*0.5*projdvgas**2
!             vissv = -alphaav*0.5*(dvdotr**2 - projepsdeltav*dvdotr + projepsdeltav*projddeltav) !-alphaav*0.5*projdvgas**2
          elseif (iav.eq.3) then
             vissv = -alphaav*0.5*projdvgasav**2
          elseif (iav.eq.2) then
             vissv = -alphaav*0.5*projdvgas**2
          else
             vissv = -alphaav*0.5*dvdotr**2
          endif
       else
          vissv = 0.
       endif
       !
       !  thermal energy terms
       !
       if (iener.eq.1) then
          vissu = 0.
       else
          vissu = alphau*(uu(i) - uu(j))
       endif
       !
       !--dissipation term in deltav
       !
       faci = 1./(dustfraci(1)*(1. - dustfraci(1)))
       facj = 1./(dustfracj(1)*(1. - dustfracj(1)))
       if (iav.eq.1 .or. iav.eq.3) then
          if (iav.eq.1) then
             !vissdv = -0.5*(projdeltavi - projdeltavj)**2
             vissdv = -0.5*dot_product(ddeltav,ddeltav)
             termdv = alphaav*vsig*rhoav1*dustfracav*grkern !*(1. - dustfracav)
          else
             vissdv = 0. ! deltav does not dissipate into thermal energy
             termdv = alphaav*vsig*rhoav1*dustfracav*grkern*(1. - dustfracav)
          endif
          if (iav.eq.3) then
             ddeltavdt(:,i) = ddeltavdt(:,i) + faci*pmassj*termdv*(-projdvgasav)*dr(:)
             ddeltavdt(:,j) = ddeltavdt(:,j) - facj*pmassi*termdv*(-projdvgasav)*dr(:)
          else
!             ddeltavdt(:,i) = ddeltavdt(:,i) + faci*pmassj*termdv*(projdeltavi - projdeltavj)*dr(:)
!             ddeltavdt(:,j) = ddeltavdt(:,j) - facj*pmassi*termdv*(projdeltavi - projdeltavj)*dr(:)
             ddeltavdt(:,i) = ddeltavdt(:,i) + faci*pmassj*termdv*(ddeltav(:))
             ddeltavdt(:,j) = ddeltavdt(:,j) - facj*pmassi*termdv*(ddeltav(:))
          endif
       elseif (iav.eq.4) then
          !
          ! this is:
          ! -rhoi/rhogasi*forcei(:)
          !
          faci = rhoi/rhogasi
          facj = rhoj/rhogasj
!          termdv = alphaav*termv
          !vsigdv = vsig !0.5*sqrt(projdeltavi**2 + projdeltavj**2)
          termdv = alphaav*vsigdv*rhoav1*grkern
          ddeltavdt(:,i) = ddeltavdt(:,i) + faci*pmassj*termdv*(projddeltav)*dr(:)
          ddeltavdt(:,j) = ddeltavdt(:,j) - facj*pmassi*termdv*(projddeltav)*dr(:)
          !vissdv = -0.5*(projdeltavi - projdeltavj)**2
          termdv = 0.
       elseif (iav.eq.2) then
          termdv = 0. ! for iav=2 the heating term comes into vissv instead of vissdv
          vissdv = 0.
          !
          !--extra dissipation term in deltav
          !
          if (projddeltav < 0.) then
             !vsigdv = 0.5*sqrt(projdeltavi**2 + projdeltavj**2)
             vsigdv = 0.5*(spsoundi + spsoundj)
             termdv = alphaav*vsigdv*rhoav1*grkern*dustfracav*(1. - dustfracav)
             ddeltavdt(:,i) = ddeltavdt(:,i) + faci*pmassj*termdv*(projddeltav)*dr(:)
             ddeltavdt(:,j) = ddeltavdt(:,j) - facj*pmassi*termdv*(projddeltav)*dr(:)
             vissdv = -0.5*projddeltav**2
          endif
       else
          termdv = 0.
          vissdv = 0.
       endif
       !
       !  add to thermal energy equation
       !
       if (damp.lt.tiny(0.)) then
          faci = rhoi/rhogasi
          facj = rhoj/rhogasj
          dudt(i) = dudt(i) + faci*pmassj*(termv*(vissv) + termu*vissu + termdv*vissdv)
          dudt(j) = dudt(j) + facj*pmassi*(termv*(vissv) - termu*vissu + termdv*vissdv)
       endif

       do k=1,ndust
          vsigeps = 0.5*(spsoundi + spsoundj)
          diffeps = alphaB*rhoav1*vsigeps*(dustfraci(k) - dustfracj(k))*grkern
          ddustevoldt(k,i) = ddustevoldt(k,i) + pmassj*diffeps
          ddustevoldt(k,j) = ddustevoldt(k,j) - pmassi*diffeps
       enddo

    endif
  end subroutine artificial_dissipation_dust

!--------------------------------------------------------------------------------------
! Artificial viscosity term for one-fluid dust in the terminal velocity approximation
! this is just regular AV term, Phantom-style, multiplied by (1 - epsilon)
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation_dust_diffusion
    implicit none
    real :: vsi,vsj,qi,qj,visc,du,cfaci,cfacj,diffu
    real :: vsigeps,alphaB,diffeps!,tstopi,tstopj

    if (dvdotr < 0.) then
       !
       ! artificial viscosity, as in PL15
       !
       vsi = alphai*spsoundi - beta*dvdotr
       vsj = alpha(1,j)*spsoundj - beta*dvdotr
       qi = -rhogasi*vsi*dvdotr
       qj = -rhogasj*vsj*dvdotr

       visc = 0.5*(qi*rho21i*grkerni + qj*rho21j*grkernj)
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
       !
       !  add to thermal energy equation
       !
       if (damp.lt.tiny(0.)) then
          !
          ! artificial thermal conductivity, as in PL15
          !
          du = uu(i) - uu(j)
          cfaci = 0.5*alphaui*rhoi*vsigu*du
          cfacj = 0.5*alpha(2,j)*rhoj*vsigu*du
          diffu = cfaci*grkerni*rho1i**2 + cfacj*grkernj*rho1j**2

          dudt(i) = dudt(i) + 0.5*qi*rho1i/rhogasi*pmassj*dvdotr*grkerni + pmassj*diffu/(1. - sum(dustfraci))
          dudt(j) = dudt(j) + 0.5*qj*rho1j/rhogasj*pmassi*dvdotr*grkernj - pmassi*diffu/(1. - sum(dustfracj))
       endif
    endif

    if (idustevol==0) then
       alphaB = 0.5*(alphaBi + alpha(3,j))
       !tstopi = get_tstop(idrag_nature,rhogasi,rhodusti,spsoundi,Kdrag)
       !tstopj = get_tstop(idrag_nature,rhogasj,rhodustj,spsoundj,Kdrag)
       !vsigeps = 0.5*(dustfraci + dustfracj)*sqrt(abs(pri - prj)*2./(rhogasi + rhogasj))
       !vsigeps = 0.5*(tstopi + tstopj)*abs(pri - prj)/rij*2./(rhogasi + rhogasj)
       !vsigeps = 4.*tstopi*tstopj/(tstopi + tstopj + 1.e-6)*abs(pri - prj)/rij*2./(rhogasi + rhogasj)
       !vsigeps = 0.5*(tstopi + tstopj)*abs(pri - prj)*2./((hi + hj)*(rhogasi + rhogasj))
       !vsigeps = 0.5*(dustfraci**2 + dustfracj**2)*0.5*(spsoundi + spsoundj)
       vsigeps = 0. !5*(dustfraci + dustfracj)*vsigu !5*(spsoundi + spsoundj)
       !vsigeps = 0.25*(dustfraci + dustfracj)*(spsoundi + spsoundj)
       do k=1,ndust
          diffeps = alphaB*rhoav1*vsigeps*(dustfraci(k) - dustfracj(k))*grkern
          ddustevoldt(k,i) = ddustevoldt(1,i) + pmassj*diffeps
          ddustevoldt(k,j) = ddustevoldt(1,j) - pmassi*diffeps
       enddo
    endif

  end subroutine artificial_dissipation_dust_diffusion

!--------------------------------------------------------------------------------------
! These are the artificial viscosity, thermal conduction terms
! Change this to change the artificial viscosity algorithm
! Inherits the values of local variables from rates
!
! This version corresponds to special relativity.
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation_sr
    implicit none
    real :: visc,alphaav,alphau
    !real :: v2i,v2j
    real :: qdiffi,qdiffj,qdiff
    real :: term,dpmomdotr,dens1j
    real :: lorentzfactorstari,lorentzfactorstarj,estari,estarj,destar
    real :: lorentzi,lorentzj,enthi,enthj,enthav,qa,qb
    real,dimension(ndimV) :: pmomstari,pmomstarj
    integer :: iavtype
    !
    !--definitions
    !
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    vsigav = max(alphaav,alphau)*vsig

    dens1i = 1./dens(i)
    dens1j = 1./dens(j)
    enthi = 1. + uu(i) + pri*dens1i
    enthj = 1. + uu(j) + prj*dens1j

    term = vsig*rhoav1*grkern

    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------

    lorentzfactorstari=1.0/(sqrt(1.0-(projvi**2)))
    lorentzfactorstarj=1.0/(sqrt(1.0-(projvj**2)))
    if (v2i.le.1.0) lorentzi=1.0/(sqrt(1.0-(v2i)))
    if (v2j.le.1.0) lorentzj=1.0/(sqrt(1.0-(v2j)))


    iavtype = 4
    select case(iavtype)
    case(4) ! qa and qb with qa*vb and qb*vb in energy equation
       ! Siegler/Riffert style
       pmomstari = veli*lorentzfactorstari !*enthi
       pmomstarj = velj*lorentzfactorstarj !*enthj
       dpmomdotr=dot_product(pmomstari-pmomstarj,dr)
       if (dpmomdotr < 0.) then
          qa = -0.5*dpmomdotr*enthi
          qb = -0.5*dpmomdotr*enthj
       else
          qa = 0.
          qb = 0.
       endif
       dpmomdotr = (qa + qb)
       destar=-alphaav*(qa*dot_product(veli,dr) + qb*dot_product(velj,dr)) &
            + alphau*(uu(i)/lorentzi - uu(j)/lorentzj)
    case(3) ! mine but with enthalpy in right place
       pmomstari = veli*lorentzfactorstari*enthi
       pmomstarj = velj*lorentzfactorstarj*enthj
       dpmomdotr=abs(dot_product(pmomstari-pmomstarj,dr))
       estari = lorentzfactorstari*enthi - pri*rho1i
       estarj = lorentzfactorstarj*enthj - prj*rho1j
       destar=alphaav*(estari - estarj)
       !destar=alphaav*(lorentzfactorstari*enthi - lorentzfactorstarj*enthj) !&
            !- alphau*(uu(i)/lorentzi - uu(j)/lorentzj)
    case(2) ! Chow/Monaghan
       pmomstari = veli*lorentzfactorstari*enthi
       pmomstarj = velj*lorentzfactorstarj*enthj
       dpmomdotr=abs(dot_product(pmomstari-pmomstarj,dr))
       estari = lorentzfactorstari*enthi - pri*rho1i
       estarj = lorentzfactorstarj*enthj - prj*rho1j
       destar=alphaav*(estari - estarj)
    case default
       ! my best one so far
       pmomstari = veli*lorentzfactorstari
       pmomstarj = velj*lorentzfactorstarj
       enthav = 0.5*(enthi + enthj)
       dpmomdotr=abs(dot_product(pmomstari-pmomstarj,dr))*enthav
       destar=alphaav*enthav*(lorentzfactorstari - lorentzfactorstarj) &
            + alphau*(uu(i)/lorentzi - uu(j)/lorentzj)
    end select

    if (dvdotr.lt.0) then
       visc = alphaav*term*dpmomdotr
       forcei(:) = (forcei(:) - pmassj*visc*dr(:))
       forcej(:) = (forcej(:) + pmassi*visc*dr(:))
    endif

    !--------------------------------------------------
    !  dissipation terms in energy equation
    !  (viscosity + resistivity + thermal conduction)
    !
    !   total energy equation
    !--------------------------------------------------
    if (iener.eq.3) then

       qdiff = 0.
       !
       !  kinetic energy terms - applied only when particles approaching
       !
       if (dvdotr.lt.0) then
          qdiff = qdiff + destar
       endif
       !
       !  add to total energy equation
       !
       dendt(i) = (dendt(i) + pmassj*term*qdiff)
       dendt(j) = (dendt(j) - pmassi*term*qdiff)

    !--------------------------------------------------
    !   thermal energy equation
    !--------------------------------------------------
    elseif (iener.eq.1) then
       qdiff = 0
       qdiffi = 0
       qdiffj = 0
       if (dvdotr.lt.0) then
          qdiffi=destar+(projvi*dpmomdotr)
          qdiffj=-destar-(projvj*dpmomdotr)
       endif

       !
       !  add to entropy equation
       !
       !if (qdiffi.gt.1.e-8 .or. qdiffj.gt.1.e-8) then
       !   print*,'error, negative entropy generation ',qdiffi,qdiffj
       !   stop
       !endif
       dudt(i) = (dudt(i) + lorentzi*pmassj*term*qdiffi)
       dudt(j) = (dudt(j) + lorentzj*pmassi*term*qdiffj)

    else
       stop 'wrong iener in special relativity'
    endif

    return
  end subroutine artificial_dissipation_sr


!--------------------------------------------------------------------------------------
! Physical viscosity terms, done using direct 2nd derivatives
! as in Espanol & Revenga (2003)
!
! SHEAR VISCOSITY ONLY
!--------------------------------------------------------------------------------------
  subroutine physical_viscosity
    implicit none
    real :: visc,grgrw

    grgrw = -2.*grkern/rij
    visc  = rhoav1*shearvisc*grgrw
    forcei(:) = forcei(:) - pmassj*visc*dvel(:)
    forcej(:) = forcej(:) + pmassi*visc*dvel(:)
    dtvisc = min(dtvisc,min(hi**2,hj**2)/shearvisc)

    if (allocated(del2v)) then
       del2v(:,i) = del2v(:,i) - pmassj*rhoav1*dvel(:)*grgrw
       del2v(:,j) = del2v(:,j) + pmassi*rhoav1*dvel(:)*grgrw
    endif

  end subroutine physical_viscosity

!----------------------------------------------------------------
! Magnetic field terms - ie force and magnetic field evolution
!
!----------------------------------------------------------------
  subroutine mhd_terms
    implicit none
    real :: dalphaterm,termi,termj
    real, dimension(ndimV) :: dBdtvisc,curlBj,curlBterm,vec,JxBxBi,JxBxBj,dBdtambi
    real :: dwdxdxi,dwdxdyi,dwdydyi,dwdxdzi,dwdydzi,dwdzdzi
    real :: dwdxdxj,dwdxdyj,dwdydyj,dwdxdzj,dwdydzj,dwdzdzj,rij1
    real :: dgradwdhi,dgradwdhj,BdotBextj,term,divBterm,etaij
    !----------------------------------------------------------------------------
    !  Lorentz force
    !----------------------------------------------------------------------------

    select case(imagforce)
    case(1)      ! vector form (dot products)

       BidotdB = dot_product(Brhoi,dB)
       BjdotdB = dot_product(Brhoj,dB)
       fmagi(:) = grkern*(dr(:)*BidotdB - projBrhoi*dB(:))
       fmagj(:) = grkern*(dr(:)*BjdotdB - projBrhoj*dB(:))

    case(3)  ! alternative symmetric formulation on aniso

       rhoij = rhoi*rhoj
       fiso = 0.5*(Brho2i*grkerni*sqrtgi + Brho2j*grkernj*sqrtgj)
       ! faniso is        B*div B        +    B dot grad B
       faniso(:) = (gdiagi(:)*Bi(:)*projBj*grkerni*sqrtgi &
                  + gdiagj(:)*Bj(:)*projBi*grkernj*sqrtgj)/rhoij
       fmagi(:) = faniso(:) - fiso*dr(:)

    case(5)      ! Morris' Hybrid form

       rhoij = rhoi*rhoj
       fiso = 0.5*(Brho2i*grkerni*sqrtgi + Brho2j*grkernj*sqrtgj)
!--note that I have tried faniso with grkerni and grkernj but much worse on mshk2
       faniso(:) = grkern*(gdiagj(:)*Bj(:)*projBj*sqrtgj &
                         - gdiagi(:)*Bi(:)*projBi*sqrtgi)/rhoij
       fmagi(:) = faniso(:) - fiso*dr(:)

    case(7)
!--vector potential force
       rij1 = 1./rij
       fiso = -1.5*(Brho2i*grkerni*sqrtgi + Brho2j*grkernj*sqrtgj)

       !--add gradh term
       fiso = fiso + (zeta(i)*rho21i*grkerni + zeta(j)*rho21j*grkernj)

       !--add B dot Bext term
       BdotBextj = dot_product(Bj(:),Bconst(:))
       fiso = fiso + 2.*(BdotBexti*rho21i*grkerni + BdotBextj*rho21j*grkernj)

       !--add dgradW/dh contribution to isotropic force term
       !  [grgrkernalti has been multiplied by 1/h^(ndim+2) already]
       dgradwdhi = -(ndim+1.)*hi1*grkernalti - rij*hi1*grgrkernalti
       dgradwdhj = -(ndim+1.)*hj1*grkernaltj - rij*hj1*grgrkernaltj
       fiso = fiso + (psi(i)*dgradwdhi + psi(j)*dgradwdhj)

       dwdxdxi = dr(1)*dr(1)*grgrkernalti + (1. - dr(1)*dr(1))*rij1*grkernalti
       if (ndim.ge.2) then
          dwdxdyi = dr(1)*dr(2)*grgrkernalti - dr(1)*dr(2)*rij1*grkernalti
          dwdydyi = dr(2)*dr(2)*grgrkernalti + (1. - dr(2)*dr(2))*rij1*grkernalti
          if (ndim.ge.3) then
             dwdxdzi = dr(1)*dr(3)*grgrkernalti - dr(1)*dr(3)*rij1*grkernalti
             dwdydzi = dr(2)*dr(3)*grgrkernalti - dr(2)*dr(3)*rij1*grkernalti
             dwdzdzi = dr(3)*dr(3)*grgrkernalti + (1. - dr(3)*dr(3))*rij1*grkernalti
          else
             dwdxdzi = 0.
             dwdydzi = 0.
             dwdzdzi = 0.
          endif
       else
          dwdxdyi = 0.
          dwdydyi = 0.
          dwdxdzi = 0.
          dwdydzi = 0.
          dwdzdzi = 0.
       endif

       dwdxdxj = dr(1)*dr(1)*grgrkernaltj + (1. - dr(1)*dr(1))*rij1*grkernaltj
       if (ndim.ge.2) then
          dwdxdyj = dr(1)*dr(2)*grgrkernaltj - dr(1)*dr(2)*rij1*grkernaltj
          dwdydyj = dr(2)*dr(2)*grgrkernaltj + (1. - dr(2)*dr(2))*rij1*grkernaltj
          if (ndim.ge.3) then
             dwdxdzj = dr(1)*dr(3)*grgrkernaltj - dr(1)*dr(3)*rij1*grkernaltj
             dwdydzj = dr(2)*dr(3)*grgrkernaltj - dr(2)*dr(3)*rij1*grkernaltj
             dwdzdzj = dr(3)*dr(3)*grgrkernaltj + (1. - dr(3)*dr(3))*rij1*grkernaltj
          else
             dwdxdzj = 0.
             dwdydzj = 0.
             dwdzdzj = 0.
          endif
       else
          dwdxdyj = 0.
          dwdydyj = 0.
          dwdxdzj = 0.
          dwdydzj = 0.
          dwdzdzj = 0.
       endif

       faniso(1) = -dBevol(3)*(Brhoi(2)*dwdxdxi*gradhi*rho1i + Brhoj(2)*dwdxdxj*gradh(j)*rho1j &
                             -(Brhoi(1)*dwdxdyi*gradhi*rho1i + Brhoj(1)*dwdxdyj*gradh(j)*rho1j)) &
                   -dBevol(2)*(Brhoi(1)*dwdxdzi*gradhi*rho1i + Brhoj(1)*dwdxdzj*gradh(j)*rho1j &
                             -(Brhoi(3)*dwdxdxi*gradhi*rho1i + Brhoj(3)*dwdxdxj*gradh(j)*rho1j)) &
                   -dBevol(1)*(Brhoi(3)*dwdxdyi*gradhi*rho1i + Brhoj(3)*dwdxdyj*gradh(j)*rho1j &
                             - Brhoi(2)*dwdxdzi*gradhi*rho1i + Brhoj(2)*dwdxdzj*gradh(j)*rho1j)

       faniso(2) = -dBevol(3)*(Brhoi(2)*dwdxdyi*gradhi*rho1i + Brhoj(2)*dwdxdyj*gradh(j)*rho1j &
                             -(Brhoi(1)*dwdydyi*gradhi*rho1i + Brhoj(1)*dwdydyj*gradh(j)*rho1j)) &
                   -dBevol(2)*(Brhoi(1)*dwdydzi*gradhi*rho1i + Brhoj(1)*dwdydzj*gradh(j)*rho1j &
                             -(Brhoi(3)*dwdxdyi*gradhi*rho1i + Brhoj(3)*dwdxdyj*gradh(j)*rho1j)) &
                   -dBevol(1)*(Brhoi(3)*dwdydyi*gradhi*rho1i + Brhoj(3)*dwdydyj*gradh(j)*rho1j &
                             - Brhoi(2)*dwdydzi*gradhi*rho1i + Brhoj(2)*dwdydzj*gradh(j)*rho1j)

       faniso(3) = -dBevol(3)*(Brhoi(2)*dwdxdzi*gradhi*rho1i + Brhoj(2)*dwdxdzj*gradh(j)*rho1j &
                             -(Brhoi(1)*dwdydzi*gradhi*rho1i + Brhoj(1)*dwdydzj*gradh(j)*rho1j)) &
                   -dBevol(2)*(Brhoi(1)*dwdzdzi*gradhi*rho1i + Brhoj(1)*dwdzdzj*gradh(j)*rho1j &
                             -(Brhoi(3)*dwdxdzi*gradhi*rho1i + Brhoj(3)*dwdxdzj*gradh(j)*rho1j)) &
                   -dBevol(1)*(Brhoi(3)*dwdydzi*gradhi*rho1i + Brhoj(3)*dwdydzj*gradh(j)*rho1j &
                             - Brhoi(2)*dwdzdzi*gradhi*rho1i + Brhoj(2)*dwdzdzj*gradh(j)*rho1j)

       !--add B Bext term
       faniso(:) = faniso(:) + (Bi(:)*rho21i*grkerni + Bj(:)*rho21j*grkernj)*projBconst

       !--add 3D term
       faniso(:) = faniso(:) - (Bevoli(:)*dot_product(curlB(:,i),dr)*rho21i*grkerni &
                                + Bevol(:,j)*dot_product(curlB(:,j),dr)*rho21j*grkernj)
       ! curlBi(:) = - (Bevoli(:)*dot_product(curlB(:,i),dr)*rho21i*grkerni &
       !                         + Bevol(:,j)*dot_product(curlB(:,j),dr)*rho21j*grkernj)
       !if (any(curlBi(:).ne.0.)) print*,i,j,curlBi(:)
       !--correct stress with constant term
       !faniso(:) = faniso(:) - (Bconst(:)*rho21i*grkerni + Bconst(:)*rho21j*grkernj)*projBconst

       !--correct with anticlumping term
       !faniso(:) = faniso(:) - stressmax*dr(:)*(rho21i*gradhi*grkernalti + rho21j*gradh(j)*grkernaltj)

       fmagi(:) = faniso(:) - fiso*dr(:)

       !--subtract B(div B) term from magnetic force
       !forcei(:) = forcei(:) - 1.0*Bi(:)*pmassj*(projBi*rho21i*grkerni + projBj*rho21j*grkernj)
       !forcej(:) = forcej(:) + 1.0*Bj(:)*pmassi*(projBi*rho21i*grkerni + projBj*rho21j*grkernj)


    case(0)
       faniso(:) = 0.
       fiso = 0.
       fmagi(:) = 0.

    case default   ! tensor formalism in generalised form
       !
       !--isotropic mag force (pressure gradient)
       !
       fiso = 0.5*(Brho2i*phij_on_phii*grkerni*sqrtgi + Brho2j*phii_on_phij*grkernj*sqrtgj)
       !
       !--anisotropic force (including variable smoothing length terms)
       !
!          faniso(:) = (Brhoi(:)*projBrhoi &
!                     - Bconst(:)*projBconst*rho21i)*phij_on_phii*grkerni &
!                    + (Brhoj(:)*projBrhoj &
!                     - Bconst(:)*projBconst*rho21j)*phii_on_phij*grkernj
       faniso(:) = (gdiagi(:)*Brhoi(:)*projBrhoi*sqrtgi &
                  - stressmax*dr(:)*rho21i)*phij_on_phii*grkerni &
                 + (gdiagj(:)*Brhoj(:)*projBrhoj*sqrtgj &
                  - stressmax*dr(:)*rho21j)*phii_on_phij*grkernj
       !
       !--add contributions to magnetic force
       !
       fmagi(:) = faniso(:) - fiso*dr(:)

    end select
    !
    !--compute rho * divergence of B
    !
    if(allocated(divBsym)) then
       divBterm = (projBi*rho21i*grkerni + projBj*rho21j*grkernj)
       divBsym(i) = divBsym(i) + pmassj*divBterm
       divBsym(j) = divBsym(j) - pmassi*divBterm
    endif
    divB(i) = divB(i) - pmassj*projdB*grkern
    divB(j) = divB(j) - pmassi*projdB*grkern
    !
    !--compute rho * current density J
    !
    dBdtambi = 0.
    if (imagforce.eq.4 .and. iambipolar.eq.0) then
    !
    !--compute J via a direct second derivative for the vector/Euler potentials
    !  (J = -\nabla^2 \alpha)
    !
       if (ndim.le.2) then
          ! this is -\nabla^2 \alpha (equivalent to -\nabla^2 A_z in 2D)
          dalphaterm = 2.*(Bevol(3,i)-Bevol(3,j))/rij
          curlB(3,i) = curlB(3,i) - pmassj*rho1j*dalphaterm*grkerni
          curlB(3,j) = curlB(3,j) + pmassi*rho1i*dalphaterm*grkernj
       else
          stop 'JxB force not yet implemented in 3D'
       endif
    else
       if (iambipolar > 0 .or. iresist==4) then
       !
       !--curl B has already been calculated, here we use it to compute
       !  the ambipolar diffusion term
       !
          if (iambipolar > 0) then
             call cross_product3D(curlB(:,i),Bi,vec) ! J x B
             call cross_product3D(vec,Bi,JxBxBi)     ! J x B x B
             call cross_product3D(curlB(:,j),Bj,vec) ! J x B
             call cross_product3D(vec,Bj,JxBxBj)     ! J x B x B
             call cross_product3D(JxBxBi,dr,curlBi)  ! misuse of curlBi variable to store (JxB) x B x \nabla
             call cross_product3D(JxBxBj,dr,curlBj)  ! misuse of curlBi variable to store (JxB) x B x \nabla
             termi = 1./(rhoi*rho_ion*gamma_ambipolar)
             termj = 1./(rhoj*rho_ion*gamma_ambipolar)
             dBdtambi(:) = termi*curlBi(:)*rho21i*grkerni + termj*curlBj(:)*rho21j*grkernj
          endif
          if (iresist==4) then
             call cross_product3D(curlB(:,i),dr,curlBi)
             call cross_product3D(curlB(:,j),dr,curlBj)
             termi = etamhd
             termj = etamhd
             dBdtambi(:) = dBdtambi(:) - termi*curlBi(:)*rho21i*grkerni - termj*curlBj(:)*rho21j*grkernj
          endif
       else
       !
       !--calculate curl B for current (only if field is 3D)
       !  this is used in the switch for the artificial resistivity term
       !
          if (imhd.gt.0 .and. iambipolar.eq.0) then ! not for vector potential or ambipolar diffusion
             call cross_product3D(dB,dr,curlBi)
             curlB(:,i) = curlB(:,i) + pmassj*curlBi(:)*grkern
             curlB(:,j) = curlB(:,j) + pmassi*curlBi(:)*grkern
          elseif (imhd.lt.0 .and. allocated(curlBsym)) then
             call cross_product3D(Bi,dr,curlBi)
             call cross_product3D(Bj(:),dr,curlBj)
             curlBterm(:) = (curlBi(:)*rho21i*grkerni + curlBj(:)*rho21j*grkernj)
             curlBsym(:,i) = curlBsym(:,i) - pmassj*curlBterm(:)
             curlBsym(:,j) = curlBsym(:,j) + pmassi*curlBterm(:)
          endif
       endif
       !
       !--add Lorentz force to total force
       !
       select case(imagforce)
       case(1)
          fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)/rho2i
          fmag(:,j) = fmag(:,j) + pmassi*fmagj(:)/rho2j
          forcei(:) = forcei(:) + pmassj*fmagi(:)/rho2i
          forcej(:) = forcej(:) + pmassi*fmagj(:)/rho2j
       case(5)    ! Morris' Hybrid force
!          fmag(:,i) = fmag(:,i) + pmassj*(faniso(:)-fiso*dr(:))
!          fmag(:,j) = fmag(:,j) + pmassi*(faniso(:)+fiso*dr(:))
          forcei(:) = forcei(:) + pmassj*(faniso(:)-fiso*dr(:))
          forcej(:) = forcej(:) + pmassi*(faniso(:)+fiso*dr(:))
       case default      ! symmetric forces fmagxi = -fmagxj
!          fmag(:,i) = fmag(:,i) + pmassj*(fmagi(:))
!          fmag(:,j) = fmag(:,j) - pmassi*(fmagi(:))
          forcei(:) = forcei(:) + pmassj*(fmagi(:))
          forcej(:) = forcej(:) - pmassi*(fmagi(:))
          if (imagforce.eq.6) then
          !--subtract B(div B) term from magnetic force
             forcei(:) = forcei(:) - Bi(:)*pmassj*(projBi*rho21i*grkerni + projBj*rho21j*grkernj)
             forcej(:) = forcej(:) + Bj(:)*pmassi*(projBi*rho21i*grkerni + projBj*rho21j*grkernj)
          endif
       end select
    endif

    !--------------------------------------------------------------------------------
    !  time derivative of magnetic field (divide by rho later)
    !  we calculate only the term for d(B/rho)/dt
    !  if evolving dB/dt the other term comes from the cty equation and
    !  is added later
    !---------------------------------------------------------------------------------
    exactderivs: if (iuse_exact_derivs.eq.0) then

    select case(imhd)
    case(4,14)   ! conservative form (explicitly symmetric)
       dBevoldti(:) = dBevoldti(:)            &
          + pmassj*(veli(:)*projBj + velj(:)*projBi)*rho1j*grkerni
       dBevoldt(:,j) = dBevoldt(:,j)          &
          - pmassi*(veli(:)*projBj + velj(:)*projBi)*rho1i*grkernj
    case(3,13)   ! goes with imagforce = 3
       dBevoldti(:) = dBevoldti(:)         &
          - pmassj*projBrhoj*dvel(:)*grkerni
       dBevoldt(:,j) = dBevoldt(:,j)       &
          - pmassi*projBrhoi*dvel(:)*grkernj
    case(2,12)   ! conservative form (no change for 1D)
       dBevoldti(:) = dBevoldti(:)            &
          - phii_on_phij*pmassj*(dvel(:)*projBrhoi - rho1i*veli(:)*projdB)*grkerni
       dBevoldt(:,j) = dBevoldt(:,j)          &
          - phij_on_phii*pmassi*(dvel(:)*projBrhoj - rho1j*velj(:)*projdB)*grkernj
    case(1,11) ! surface flux-conservative (usual) form
       dBevoldti(:) = dBevoldti(:)            &
          - phii_on_phij*pmassj*((dvel(:)*projBrhoi)*grkerni + rhoi*dBdtambi(:))
       dBevoldt(:,j) = dBevoldt(:,j)          &
          - phij_on_phii*pmassi*((dvel(:)*projBrhoj)*grkernj - rhoj*dBdtambi(:))
    case(-1) ! vector potential evolution - crap gauge
       termi = pmassj*projvi*grkerni
       dBevoldti(:) = dBevoldti(:) - termi*dBevol(:)
       termj = pmassi*projvj*grkernj
       dBevoldt(:,j) = dBevoldt(:,j) - termj*dBevol(:)
    case(-2) ! vector potential evolution - Axel gauge
       termi = pmassj*dot_product(Bevol(:,i),dvel(:))*grkerni
       dBevoldti(:) = dBevoldti(:) + termi*dr(:)
       termj = pmassi*dot_product(Bevol(:,j),dvel(:))*grkernj
       dBevoldt(:,j) = dBevoldt(:,j) + termj*dr(:)
    end select

    endif exactderivs

    !--------------------------------------------
    !  real resistivity in induction equation
    !--------------------------------------------
    if (iresist.gt.0 .and. (iresist.ne.2 .and. iresist.ne.4)) then
       etaij = 0.5*(etai + etaj)
       !etaij = 2.*etai*etaj/(etai + etaj)
       if (imhd.gt.0) then
          dBdtvisc(:) = -2.*etaij*dB(:)/(rij + epsilon(rij))
       else !--vector potential resistivity
          dBdtvisc(:) = 2.*etaij*dBevol(:)/rij
       endif
       !
       !--add to dB/dt (converted to d(B/rho)/dt later if required)
       !
       dBevoldti(:) = dBevoldti(:) - rhoi*pmassj*0.5*(rho1i**2*grkerni + rho1j**2*grkernj)*dBdtvisc(:)
       dBevoldt(:,j) = dBevoldt(:,j) + rhoj*pmassi*0.5*(rho1i**2*grkerni + rho1j**2*grkernj)*dBdtvisc(:)
!       dBevoldti(:) = dBevoldti(:) - pmassj*rho1j*dBdtvisc(:)*grkerni
!       dBevoldt(:,j) = dBevoldt(:,j) + pmassi*rho1i*dBdtvisc(:)*grkernj
       !
       !  add to thermal energy equation
       !
       if (iener.eq.3 .or. iener.eq.1) then
          stop 'unimplemented iener for physical resistivity'
       elseif (iener.gt.0) then
          term = -etaij*rho1i*rho1j*dot_product(dB,dB)*grkern/rij
          dudt(i) = dudt(i) + pmassj*term
          dudt(j) = dudt(j) + pmassi*term
       endif
    endif

    if (idivBzero.ge.2) then ! add hyperbolic correction term
       !gradpsiterm = (psi(i)-psi(j))*grkern ! (-ve grad psi)
       gradpsiterm = (psi(i)*rho21i*grkerni + psi(j)*rho21j*grkernj)
       gradpsi(:,i) = gradpsi(:,i) - pmassj*gradpsiterm*dr(:)
       gradpsi(:,j) = gradpsi(:,j) + pmassi*gradpsiterm*dr(:)
    endif

    return
  end subroutine mhd_terms

!----------------------------------------------------------------
!  Derivatives of dust to gas ratio and deltav
!  for (full) one fluid dust, as in Laibe & Price (2014b)
!----------------------------------------------------------------
  subroutine dust_derivs
    implicit none
    real :: termi,termj,term,dterm,du,si,sj,deps
    real, dimension(ndimV) :: prdustterm,dtermvec,termivec,termjvec

    !
    !--time derivative of dust to gas ratio
    !  (symmetric derivative to conserve dust/gas mass)
    !
    do k=1,ndust
       select case(idustevol)
       case(2)
          si = sqrt(rhogasi*rhodusti(k))
          sj = sqrt(rhogasj*rhodustj(k))
          termi = (1. - dustfraci(k))*rho1i*projdeltavi*grkerni
          termj = (1. - dustfracj(k))*rho1j*projdeltavj*grkernj
          term = termi + termj
          ddustevoldt(k,i) = ddustevoldt(k,i) - pmassj*sj*term
          ddustevoldt(k,j) = ddustevoldt(k,j) + pmassi*si*term
       case(1)
          si = sqrt(rhodusti(k))
          sj = sqrt(rhodustj(k))
          termi = (1. - dustfraci(k))*rho1i*projdeltavi*grkerni
          termj = (1. - dustfracj(k))*rho1j*projdeltavj*grkernj
          term = termi + termj
          ddustevoldt(k,i) = ddustevoldt(k,i) - pmassj*sj*term
          ddustevoldt(k,j) = ddustevoldt(k,j) + pmassi*si*term
       case default
          termi = rhogrhodonrhoi(k)*projdeltavi*rho21i*grkerni
          termj = rhogrhodonrhoj(k)*projdeltavj*rho21j*grkernj
          term  = termi + termj
          ddustevoldt(k,i) = ddustevoldt(k,i) - pmassj*term
          ddustevoldt(k,j) = ddustevoldt(k,j) + pmassi*term
       end select
    enddo

    !
    !--time derivative of deltav
    !  (here only bits that involve sums over particles, i.e. not decay term)
    !
    !--high mach number term
    termi    = (rhogasi - rhodusti(1))*rho1i*deltav2i
    termj    = (rhogasj - rhodustj(1))*rho1j*deltav2j
    dterm    = 0.5*(termi - termj)
    !
    !--note: need to add term to d/dt(deltav) using the gas-only force
    !  i.e. before we add the anisotropic pressure term to the forces.
    !  Using forcei and forcej directly means that we automatically include viscosity, MHD etc.
    !
    ! ORIGINAL LP14b VERSION
    ddeltavdt(:,i) = ddeltavdt(:,i) + rho1i*pmassj*(dvel(:)*projdeltavi + dterm*dr(:))*grkerni ! add forcei term once finished
    ddeltavdt(:,j) = ddeltavdt(:,j) + rho1j*pmassi*(dvel(:)*projdeltavj + dterm*dr(:))*grkernj - rhoj/rhogasj*forcej(:)

    ! TIMOTHEE VERSION
    !termivec = (rhogasi - rhodusti(1))*rho1i*deltavi(:)
    !termjvec = (rhogasj - rhodustj(1))*rho1j*deltavj(:)
    !dtermvec = (termivec - termjvec)
    !deps = dustfraci(1) - dustfracj(1)
    !ddeltavdt(:,i) = ddeltavdt(:,i) + rho1i*pmassj*(dvel(:)*projdeltavi + dtermvec(:)*projdeltavi &
    !               + deltavi(:)*deps*projdeltavi)*grkerni ! add forcei term once finished
    !ddeltavdt(:,j) = ddeltavdt(:,j) + rho1j*pmassi*(dvel(:)*projdeltavj + dtermvec(:)*projdeltavj &
    !               + deltavj(:)*deps*projdeltavj)*grkernj - rhoj/rhogasj*forcej(:)

    !
    !--anisotropic pressure term
    !
    prdustterm(:) = rhogrhodonrhoi(1)*deltavi(:)*projdeltavi*rho21i*grkerni &
                  + rhogrhodonrhoj(1)*deltavj(:)*projdeltavj*rho21j*grkernj

    fextrai(:) = fextrai(:) - pmassj*(prdustterm(:))
    fextraj(:) = fextraj(:) + pmassi*(prdustterm(:))

    !
    !--thermal energy equation: add Pg/rhog*div(vgas) and deltav.grad(u) term
    !
    if (iener.gt.0) then
       du = uu(i) - uu(j)
       dudt(i) = dudt(i) + pmassj*(pri*rho1i/rhogasi*projdvgas - rhodusti(1)*rho21i*du*projdeltavi)*grkerni
       dudt(j) = dudt(j) + pmassi*(prj*rho1j/rhogasj*projdvgas - rhodustj(1)*rho21j*du*projdeltavj)*grkernj
    endif

  end subroutine dust_derivs

!----------------------------------------------------------------
!  Derivative of dust fraction and thermal energy
!  for one fluid dust in the terminal velocity approximation
!----------------------------------------------------------------
  subroutine dust_derivs_diffusion
    real :: diffterm, Di, Dj, du, Dav
    real :: tstopi, tstopj, pdvtermi, pdvtermj
    real :: si, sj, grgrkern
    integer :: k

    do k=1,ndust
       tstopi = get_tstop(idrag_nature,rhogasi,rhodusti(k),spsoundi,Kdrag)
       tstopj = get_tstop(idrag_nature,rhogasj,rhodustj(k),spsoundj,Kdrag)

       select case(idustevol)
       case(2,4)
          Di = rho1i*tstopi/rhogasi
          Dj = rho1j*tstopj/rhogasj
          Dav = 0.5*(Di + Dj)
          si = sqrt(rhogasi*rhodusti(k))
          !print*,' si = ',si,dustevol(i)
          sj = sqrt(rhogasj*rhodustj(k))
          grgrkern = -2.*grkern/rij
          diffterm = Dav*(pri - prj)*grgrkern
          ddustevoldt(k,i) = ddustevoldt(k,i) + pmassj*sj*rho1j*diffterm
          ddustevoldt(k,j) = ddustevoldt(k,j) - pmassi*si*rho1i*diffterm
       case(1) ! sqrt(rho*eps)
          Di = rho1i*tstopi
          Dj = rho1j*tstopj
          Dav = 0.5*(Di + Dj)
          si = sqrt(rhodusti(k))
          !print*,' si = ',si,dustevol(i)
          sj = sqrt(rhodustj(k))
          grgrkern = -2.*grkern/rij
          diffterm = rho1i*rho1j*Dav*(pri - prj)*grgrkern
          ddustevoldt(k,i) = ddustevoldt(k,i) + pmassj*sj*diffterm
          ddustevoldt(k,j) = ddustevoldt(k,j) - pmassi*si*diffterm
       case default
          Di = rhodusti(k)*rho1i*tstopi
          Dj = rhodustj(k)*rho1j*tstopj
          !Di = dustfraci*tstopi
          !Dj = dustfracj*tstopj
          !if (Di + Dj > 0.) then
             Dav = 0.5*(Di + Dj)
             !Dav = 2.*Di*Dj/(Di + Dj)
             !Dav = sqrt(0.5*(Di**2 + Dj**2)) !0.5*(Di + Dj) !2.*Di*Dj/(Di + Dj)
          !else
          !   Dav = 0.
          !endif
          grgrkern = -2.*grkern/rij
          !grgrkern = 0.5*(grgrkernalti + grgrkernaltj)

          diffterm = rho1i*rho1j*Dav*(pri - prj)*grgrkern
          ddustevoldt(k,i) = ddustevoldt(k,i) + pmassj*diffterm
          ddustevoldt(k,j) = ddustevoldt(k,j) - pmassi*diffterm
       end select
    enddo

    !
    !--thermal energy equation: add Pg/rhog*div(vgas) and deltav.grad(u) term
    !
    if (iener.gt.0) then
       du = uu(i) - uu(j)
       pdvtermi = pri*rho1i/rhogasi*pmassj*dvdotr*grkerni
       pdvtermj = prj*rho1j/rhogasj*pmassi*dvdotr*grkernj
       dudt(i) = dudt(i) + pdvtermi - 0.5*pmassj*du*diffterm/(1. - sum(dustfraci))
       dudt(j) = dudt(j) + pdvtermj - 0.5*pmassi*du*diffterm/(1. - sum(dustfracj))
    endif

  end subroutine dust_derivs_diffusion

!---------------------------------------------------------------------------
!  quantum terms
!
!---------------------------------------------------------------------------
  subroutine quantum_terms
    real :: qsphterm(ndim),termi(ndim),termj(ndim)

    do k=1,ndim
       termi(k) = dot_product(P_Q(k,:,i),dr(1:ndim))
       termj(k) = dot_product(P_Q(k,:,j),dr(1:ndim))
    enddo

    qsphterm(:) = termi(:)*rho21i*grkerni + termj(:)*rho21j*grkernj

    forcei(1:ndim) = forcei(1:ndim) - pmassj*qsphterm(:)
    forcej(1:ndim) = forcej(1:ndim) + pmassi*qsphterm(:)

  end subroutine quantum_terms

!----------------------------------------------------------------------------
!  energy equation if evolving the total energy/particle
!
!----------------------------------------------------------------------------
  subroutine energy_equation
    implicit none
    real, dimension(ndimV) :: prvterm
    real :: Brhoidotvj,Brhojdotvi,Brhoidotdv, Brhojdotdv
    real :: prvaniso, projprv, prvanisoi, prvanisoj

    Brhoidotvj = dot_product(Brhoi,velj)
    Brhojdotvi = dot_product(Brhoj,veli)
    !
    ! (isotropic stress)
    !
    prvterm(:) = (Prho2i+0.5*Brho2i)*phii_on_phij*velj(:)*sqrtgi*grkerni &
               + (Prho2j+0.5*Brho2j)*phij_on_phii*veli(:)*sqrtgj*grkernj
    projprv = dot_product(prvterm,dr)
    !
    ! (anisotropic stress)
    !
    if (imhd.ne.0 .and. imagforce.eq.2) then
       !
       !  with anticlumping term must add contribution from anisotropic term
       !  explicitly and subtract contribution from B/rho version of induction equation
       !  (some of these terms cancel if the kernels are the same)
       !
       Brhoidotdv = dot_product(Brhoi,dvel) ! these terms are from the B/rho induction eq.
       Brhojdotdv = -dot_product(Brhoj,dvel) ! -ve so velj-veli instead of veli-velj

       prvanisoi = -(dot_product(veli,faniso) - Brhoidotdv*projBrhoi*phii_on_phij*grkerni)
       prvanisoj = -(dot_product(velj,faniso) - Brhojdotdv*projBrhoj*phij_on_phii*grkernj)

    else
       !
       !  without anticlumping term
       !
       prvaniso =  - Brhoidotvj*projBrhoi*phii_on_phij*grkerni      &
                   - Brhojdotvi*projBrhoj*phij_on_phii*grkernj
       prvanisoi = prvaniso
       prvanisoj = prvaniso
    endif

    dendt(i) = dendt(i) - pmassj*(projprv+prvanisoi)
    dendt(j) = dendt(j) + pmassi*(projprv+prvanisoj)
    !
    ! (source term for hyperbolic divergence correction)
    !
    if (idivBzero.ge.2) then
       dendt(i) = dendt(i) + pmassj*projBi*gradpsiterm
       dendt(j) = dendt(j) + pmassi*projBj*gradpsiterm
    endif

    return
  end subroutine energy_equation

end subroutine get_rates
