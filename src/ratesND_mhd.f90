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

!!--------------------------------------------------------------------
!! Computes the rates of change of the conserved variables
!! (forces, energy etc)
!! This is the core of the SPH algorithm
!!--------------------------------------------------------------------

subroutine get_rates
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
 use utils, only:cross_product3D
  
 use fmagarray
 use derivB
 use get_neighbour_lists
 use externf,       only:external_forces,pequil
 use dust,          only:get_tstop
!
!--define local variables
!
 implicit none
 integer :: i,j,n
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
 integer :: itypei
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
 real :: BidotdB,BjdotdB,Brho2i,Brho2j,BdotBexti,etai
 real :: projBrhoi,projBrhoj,projBi,projBj,projdB,projBconst
 real, dimension(:,:), allocatable :: curlBsym
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
 real, dimension(:,:), allocatable :: graddivv
!
!  (alternative forms)
!
 real, dimension(:), allocatable :: phi
 real :: phii,phii1,phii_on_phij,phij_on_phii,uui,uuj
!
!  (one fluid dust)
!
 real, dimension(ndimV) :: deltavi,deltavj
 real :: dustfraci,dustfracj
 real :: rhogrhodonrhoi, rhogrhodonrhoj
 real :: deltav2i,deltav2j
 real :: dtstop,projdeltavi,projdeltavj
 real :: rhogasi,rhodusti,rhogasj,rhodustj,projdvgas
 real, dimension(ndimV) :: vgasi,vgasj,dvgas,fextrai,fextraj
 real :: sum,tstop,ratio,dvmax

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
 allocate ( phi(ntotal), del2u(ntotal), graddivv(ndimV,ntotal), STAT=ierr )
 if (ierr.ne.0) write(iprint,*) ' Error allocating phi, ierr = ',ierr
 if (imhd.lt.0) then
    allocate( curlBsym(ndimV,ntotal), STAT=ierr )
    if (ierr.ne.0) write(iprint,*) ' Error allocating curlBsym, ierr = ',ierr
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
  if (imhd.gt.0 .and. iambipolar.eq.0) curlB(:,i) = 0.0
  if (allocated(curlBsym)) curlBsym(:,i) = 0.
  if (allocated(del2v)) del2v(i) = 0.
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
  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     ddustevoldt(i) = 0.
     ddeltavdt(:,i) = 0.
  endif
 enddo

!
!--calculate maximum neg stress for instability correction
!  
 stressmax = 0.
 if (imhd.ne.0 .and. (imagforce.eq.2 .or. imagforce.eq.7)) then
    do i=1,ntotal
       B2i = dot_product(Bfield(:,i),Bfield(:,i))
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
       prneti = pri - pequil(iexternal_force,xi(:),rhoi)
       pmassi = pmass(i)
       Prho2i = pri*rho21i
       spsoundi = spsound(i)
       uui = uu(i)
       veli(:) = vel(:,i)
       v2i = dot_product(veli(:),veli(:))
       alphai = alpha(1,i)
       alphaui = alpha(2,i)
       phii = phi(i)
       phii1 = 1./phii
       sqrtgi = sqrtg(i)
       gdiagi(:) = 1.

       ! one fluid dust definitions
       if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
          dustfraci  = dustfrac(i)
          rhodusti       = dustfraci*rhoi          
          rhogasi        = rhoi - rhodusti
          deltavi(:)     = deltav(:,i)
          deltav2i       = dot_product(deltavi,deltavi)
          rhogrhodonrhoi = rhogasi*rhodusti*rho1i
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
          B2i = dot_product(Bi,Bi)
          Brho2i = B2i*rho21i
          valfven2i = B2i*rho1i
          alphaBi = alpha(3,i)
          if (imhd.lt.0) Bevoli(:) = Bevol(:,i)
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
              if (rij.le.1.e-8 .and. itype(j).eq.itypei) then
                 nclumped = nclumped + 1
                 if (rij.lt.tiny(rij)) then
                    write(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                    call quit
                 endif
              elseif (rij.le.epsilon(rij)) then
                 dr = 0.  ! can happen with gas/dust if particles on top of each other
              else
                 dr(1:ndim) = dx(1:ndim)/rij      ! unit vector
              endif
              if (itype(j).eq.itypei &
                  .or.(itype(j).eq.itypebnd .or. itype(j).eq.itypebnd2) &
                  .or.(itypei  .eq.itypebnd .or. itype(j).eq.itypebnd2)) then              
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
 endif
 if (itiming) call cpu_time(t4)

 if (trace) write(iprint,*) 'Finished main rates loop'
 fhmax = 0.0
 dtforce = 1.e6
!! dtcourant2 = 1.e6

!----------------------------------------------------------------------------
!  loop over the particles again, subtracting external forces and source terms
!----------------------------------------------------------------------------
 
!
!--calculate maximum vsig over all the particles for use in the hyperbolic cleaning
!
 if (imhd.ne.0 .and. idivBzero.ge.2) then
    vsig2max = vsigmax**2
 endif
 
 fmean(:) = 0.
 dtdrag = huge(dtdrag)
 sum = 0.
 ratio = 0.
 if (idust.eq.2 .and. h_on_csts_max > 1.) then
    print*,' WARNING: violating h < cs*ts resolution criterion by factor of ',h_on_csts_max
 endif
 do i=1,npart
 
    rhoi  = rho(i)
    rho1i = 1./rhoi
    if (ivisc.gt.0 .and. idust.ne.1 .and. idust.ne.4) then
       sum = sum + pmass(i)*dot_product(vel(:,i),force(:,i))
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
       dustfraci = dustfrac(i)
       rhodusti = rhoi*dustfraci       
       rhogasi  = rhoi - rhodusti
       
       tstop = get_tstop(idrag_nature,rhogasi,rhodusti,spsound(i),Kdrag)
       dtdrag = min(dtdrag,tstop)
       !
       !--d/dt(deltav)  : add terms that do not involve sums over particles
       !
       if (dustfraci.gt.0.) then
          dtstop   = 1./tstop
          ddeltavdt(:,i) = ddeltavdt(:,i) - deltav(:,i)*dtstop
       else
          dtstop = 0.
          ddeltavdt(:,i) = 0.
       endif
       
       !
       !--du/dt: add thermal energy dissipation from drag
       !  (maybe should do this in the step routine along with the implicit drag in ddeltav/dt?)
       !
       if (iener.gt.0) then
          deltav2i = dot_product(deltav(:,i),deltav(:,i))
          dudt(i) = dudt(i) + rhodusti*rho1i*deltav2i*dtstop
       endif
    elseif (idust.eq.3 .or. idust.eq.4) then
       !--------------------------------------------
       !  one fluid dust in diffusion approximation
       !--------------------------------------------
       dustfraci = dustfrac(i)
       rhodusti = rhoi*dustfraci
       rhogasi  = rhoi - rhodusti
       !
       !--compute stopping time for drag timestep
       !
       tstop = get_tstop(idrag_nature,rhogasi,rhodusti,spsound(i),Kdrag)
       ! CAUTION: Line below must be done BEFORE external forces have been applied
       deltav(:,i) = -1./(1. - dustfraci)*force(:,i)*tstop
       ratio = max(dustfraci*tstop/dtcourant,ratio)
       dtdrag = min(dtdrag,0.5*hh(i)**2/(dustfraci*tstop*spsound(i)**2))
       dvmax = maxval(abs(deltav(:,i)))
       if (dvmax > 0.) then
          dtdrag = min(dtdrag,0.1*hh(i)/dvmax)
       endif
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
          elseif (imhd.gt.0) then ! not for vector potential
             curlB(:,i) = curlB(:,i)*rho1i
          endif
       endif
       divB(i) = divB(i)*rho1i
    endif
!
!--add external (body) forces
!
    if (iexternal_force.ne.0) then
       call external_forces(iexternal_force,x(1:ndim,i),fexternal(1:ndimV), &
                            ndim,ndimV,vel(1:ndimV,i),hh(i),spsound(i))
       force(1:ndimV,i) = force(1:ndimV,i) + fexternal(1:ndimV)
    endif
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
!--calculate resistive timestep (bootstrap onto force timestep)
!
    if (iresist.gt.0 .and. etamhd.gt.tiny(etamhd)) then
       fhmax = max(fhmax,etamhd/hh(i)**2)
    endif

!
!--calculate simpler estimate of vsig for divergence cleaning and 
!  in the dissipation switches
!
    valfven2i = 0.
    if (imhd.ne.0) valfven2i = dot_product(Bfield(:,i),Bfield(:,i))/dens(i)       
    vsig2 = spsound(i)**2 + valfven2i     ! approximate vsig only
    vsig = SQRT(vsig2)
!
!--if evolving B instead of B/rho, add the extra term from the continuity eqn
!  (note that the dBevoldt term should be divided by rho)
!  
    if (imhd.ge.11) then  ! evolving B
       dBevoldt(:,i) = sqrtg(i)*dBevoldt(:,i) + Bevol(:,i)*rho1i*drhodt(i)
       if (idivBzero.ge.2) then
          gradpsi(:,i) = gradpsi(:,i)*rho1i
       endif
    elseif (imhd.gt.0) then ! evolving B/rho
       dBevoldt(:,i) = sqrtg(i)*dBevoldt(:,i)*rho1i
       if (idivBzero.ge.2) then
          gradpsi(:,i) = gradpsi(:,i)*rho1i**2
       endif
    elseif (imhd.eq.-1) then ! vector potential evolution, crap gauge
       !
       !--add the v x B term
       !
       call cross_product3D(vel(:,i),Bfield(:,i),curlBi)
       dBevoldt(:,i) = dBevoldt(:,i)*rho1i + curlBi(:)
    elseif (imhd.eq.-2) then ! vector potential evolution, Axel gauge
       !
       !--get v x Bext
       !
       call cross_product3D(vel(:,i),Bconst(:),curlBi)
       !--add v x Bext plus the existing term, which includes dissipation
       dBevoldt(:,i) = curlBi(:) + dBevoldt(:,i)*rho1i
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
    else
       dBevoldt(:,i) = 0.
    endif
!
!--if using the thermal energy equation, set the energy derivative
!  (note that dissipative terms are calculated in rates, but otherwise comes straight from cty)
!

    if (iener.eq.3) then
       ! could do this in principle but does not work with
       ! faniso modified by subtraction of Bconst
       !dendt(i) = dot_product(vel(:,i),force(:,i)) + dudt(i) &
       !         + 0.5*(dot_product(Bfield(:,i),Bfield(:,i))*rho1i**2) &
       !         + dot_product(Bfield(:,i),dBevoldt(:,i))*rho1i
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
       sum = sum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - dot_product(vel(:,i),fexternal(:)) &
             + rhogasi*rhodusti*rho1i**2*dot_product(deltav(:,i),ddeltavdt(:,i)) &
             + ((1. - 2.*dustfraci)*0.5*dot_product(deltav(:,i),deltav(:,i)) - uu(i))*ddustevoldt(i) &
             + rhogasi*rho1i*dudt(i))
    elseif (idust.eq.4) then
       sum = sum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
             - dot_product(vel(:,i),fexternal(:)) &
             - uu(i)*ddustevoldt(i) &
             + (1. - dustfrac(i))*dudt(i))
    endif
 enddo
 if ((idust.eq.1 .or. idust.eq.4) .and. abs(sum).gt.epsilon(sum) .and. (iener.ge.2)) then
    !print*,' SUM (should be zero if conserving energy) = ',sum
    !read*
 endif
 if (idust.ne.0 .and. ratio > 1.) print "(a,g8.3,a)",' WARNING: max ts/dt = ',ratio,' approximation not valid'

 if (ivisc.gt.0) print*,' dEk/dt = ',sum
 
 if (imhd.ne.0 .and. &
    (sqrt(dot_product(fmean,fmean)).gt.1.e-8) &
     .and. mod(nsteps,100).eq.0) print*,'WARNING: fmean = ',fmean(:)
!
!--calculate timestep constraint from the forces
!  dtforce is returned together with dtcourant to the main timestepping loop
!
 if (fhmax.lt.0.) then
    write(iprint,*) 'rates: fhmax <=0 :',fhmax
    call quit
 elseif (fhmax.gt.0.) then
    if (dtforce.gt.0.0) dtforce = sqrt(1./fhmax)    
 else
    dtforce = 1.e6
 endif    
 !!print*,'dtcourant = ',dtcourant,dtcourant2,0.2*dtcourant2
 !!dtcourant = 0.2*dtcourant2
!
!--set rates to zero on ghosts/fixed particles
!
 do i=1,ntotal      ! using ntotal just makes sure they are zero for ghosts
    if (itype(i).eq.itypebnd .or. i.gt.npart) then
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
 if (allocated(phi)) deallocate(phi,del2u,graddivv)
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

!----------------------------------------------
! Evaluate drag forces on a pair of particles
!----------------------------------------------
  subroutine drag_forces
    use kernels, only:interpolate_kerneldrag,interpolate_kernel  
    use options, only:idrag_nature,Kdrag
    use dust,    only:get_tstop
    implicit none
    integer :: itypej
    real    :: dv2,vij,projv,dragterm,dragterm_en
    real    :: spsoundgas,ts
    real    :: wabj,wab,hfacwabj
    real, dimension(ndimV) :: drdrag
    real, parameter :: pi  = 3.141592653589
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
    hfacwabj = (1./hh(j)**ndim)
    call interpolate_kerneldrag(q2i,wabi)
    wabi     = wabi*hfacwabi
    call interpolate_kerneldrag(q2j,wabj)
    wabj     = wabj*hfacwabj
!
!--get particle j properties
!
    itypej = itype(j)
    pmassj = pmass(j)
    rhoj   = rho(j)
!
!--calculate the stopping time for use in the drag force
!
    projv = dot_product(dvel,drdrag)
    if (itypei.eq.itypegas) then
       spsoundgas = spsound(i)
       wab = wabi
       ts = get_tstop(idrag_nature,rhoi,rhoj,spsoundgas,Kdrag)
       h_on_csts_max = max(h_on_csts_max,hh(i)/(spsoundgas*ts))
    else
       if (itypej.ne.itypegas) return
       spsoundgas = spsound(j)
       wab = wabj
       ts = get_tstop(idrag_nature,rhoj,rhoi,spsoundgas,Kdrag)
       h_on_csts_max = max(h_on_csts_max,hh(j)/(spsoundgas*ts))
    endif
    ts_min = min(ts_min,ts)
!
!--update the force and the energy
!
    dragterm    = ndim*wab/((rhoi + rhoj)*ts)*projv
    dragterm_en = dragterm*projv
    forcei(:)   = forcei(:)  - dragterm*pmassj*drdrag(:)   
    force(:,j)  = force(:,j) + dragterm*pmassi*drdrag(:)
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
    implicit none
    real :: prstar,prstari,prstarj,vstar
    real :: dvsigdtc
    real :: hav,hav1,h21,q2
    real :: hfacwab,hfacwabj,hfacgrkern,hfacgrkernj
    real :: wabalti,wabaltj,wabalt
    real :: altrhoi,altrhoj
    real :: term

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
    projvi = dot_product(veli,dr)
    projvj = dot_product(velj,dr)
    rhoj = rho(j)
    rho1j = 1./rhoj
    rho2j = rhoj*rhoj
    rho21j = rho1j*rho1j
!!    rhoj5 = sqrt(rhoj)
    rhoij = rhoi*rhoj
    rhoav1 = 0.5*(rho1i + rho1j)   !2./(rhoi + rhoj)
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       dustfracj  = dustfrac(j)
       rhodustj       = rhoj*dustfracj
       rhogasj        = rhoj - rhodustj
       rhogrhodonrhoj = rhogasj*rhodustj*rho1j
       deltavj(:)     = deltav(:,j)
       deltav2j       = dot_product(deltavj,deltavj)
       vgasi(:)       = veli(:) - dustfraci*deltavi(:)
       vgasj(:)       = velj(:) - dustfracj*deltavj(:)
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
    gdiagj(:) = 1.
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
       B2j = dot_product(Bj,Bj)
       valfven2j = B2j*rho1j
       Brho2j = B2j*rho21j
       projBconst = dot_product(Bconst,dr)
       if (imhd.lt.0) then
          dBevol(:) = Bevoli(:) - Bevol(:,j)
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
       vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
       vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))
       vsigB = 0.5*(vsigi + vsigj) + abs(dvdotr)
       !vsigB = sqrt(dot_product(dvel - dvdotr,dvel - dvdotr))
       !vsigB = 0.5*(sqrt(valfven2i) + sqrt(valfven2j))
    else
       vsigi = spsoundi
       vsigj = spsoundj
       vsigB = 0.
    endif

    vsig = 0.5*(max(vsigi + vsigj - beta*dvdotr,0.0)) ! also used where dvdotr>0 in MHD
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
       if (idust.eq.1) then
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
          prterm = prstari*phii_on_phij*rho21i*sqrtgi*grkerni &
                 + prstarj*phij_on_phii*rho21j*sqrtgj*grkernj
       else
          prterm = phii_on_phij*Prho2i*sqrtgi*grkerni &
                 + phij_on_phii*Prho2j*sqrtgj*grkernj
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

!------------------------------------------------------------------------
!  total energy equation (thermal energy equation terms calculated
!                         outside loop and in artificial_dissipation)
!------------------------------------------------------------------------
   
    if (iener.eq.3) then
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
    endif
!
!  Continuity equation if density not done by summation
!
    if (icty.ge.1) then
       drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni
       drhodt(j) = drhodt(j) + pmassi*dvdotr*grkernj
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
    dpmomdotr = -dvdotr
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
    real :: dudti,dudtj

    dudti = 0.
    dudtj = 0.
    if (dvdotr < 0.) then
       !
       ! artificial viscosity, as in Phantom
       !
       vsi = max(alphai*spsoundi - beta*dvdotr,0.)
       vsj = max(alpha(1,j)*spsoundj - beta*dvdotr,0.)
       qi = -rhoi*vsi*dvdotr
       qj = -rhoj*vsj*dvdotr

       visc = 0.5*(qi*rho21i*grkerni + qj*rho21j*grkernj)
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
       !
       !  add to thermal energy equation
       !
       if (damp.lt.tiny(0.)) then
          dudti = 0.5*qi*rho21i*pmassj*dvdotr*grkerni
          dudtj = 0.5*qj*rho21j*pmassi*dvdotr*grkernj
       endif
    endif

    if (damp < tiny(damp)) then
       !
       ! artificial thermal conductivity, as in PL15
       !
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
    real :: visc,alphaav,alphau
    real :: vissv,vissu,vissdv,vsigdv
    real :: term,dpmomdotr,faci,facj
    real :: termu,termv,termdv
    
    real :: dustfracav,projdvgasav,projddeltav,projepsdeltav
    real, dimension(ndimV) :: ddeltav
    logical :: allpairs
    !
    !--definitions
    !
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    vsigav = max(alphaav,alphau)*vsig

    dustfracav = 0.5*(dustfraci + dustfracj)
    projdvgasav = dvdotr - dustfracav*(projdeltavi - projdeltavj)
    ddeltav     = deltavi(:) - deltavj(:)
    projddeltav = dot_product(ddeltav,dr)
    projepsdeltav = dustfraci*projdeltavi - dustfracj*projdeltavj

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
       faci = 1./(dustfraci*(1. - dustfraci))
       facj = 1./(dustfracj*(1. - dustfracj))
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
       
    endif
  end subroutine artificial_dissipation_dust

!--------------------------------------------------------------------------------------
! Artificial viscosity term for one-fluid dust in the terminal velocity approximation
! this is just regular AV term, Phantom-style, multiplied by (1 - epsilon)
!--------------------------------------------------------------------------------------
  subroutine artificial_dissipation_dust_diffusion
    implicit none
    real :: vsi,vsj,qi,qj,visc,du,cfaci,cfacj,diffu
    !real :: tstopi,tstopj,vsigeps,alphaB,diffeps

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
          
          dudt(i) = dudt(i) + 0.5*qi*rho1i/rhogasi*pmassj*dvdotr*grkerni + pmassj*diffu/(1. - dustfraci)
          dudt(j) = dudt(j) + 0.5*qj*rho1j/rhogasj*pmassi*dvdotr*grkernj - pmassi*diffu/(1. - dustfracj)
       endif
    endif
    
    if (.not.use_sqrtdustfrac) then
       !alphaB = 0.5*(alphaBi + alpha(3,j))
       !tstopi = get_tstop(idrag_nature,rhogasi,rhodusti,spsoundi,Kdrag)
       !tstopj = get_tstop(idrag_nature,rhogasj,rhodustj,spsoundj,Kdrag)
       !vsigeps = 0.5*(dustfraci + dustfracj)*sqrt(abs(pri - prj)*2./(rhogasi + rhogasj))
       !vsigeps = 0.5*(tstopi + tstopj)*abs(pri - prj)/rij*2./(rhogasi + rhogasj)
       !vsigeps = 4.*tstopi*tstopj/(tstopi + tstopj + 1.e-6)*abs(pri - prj)/rij*2./(rhogasi + rhogasj)
       !vsigeps = 0.5*(tstopi + tstopj)*abs(pri - prj)*2./((hi + hj)*(rhogasi + rhogasj))
       !vsigeps = 0.5*(dustfraci**2 + dustfracj**2)*0.5*(spsoundi + spsoundj)
       !vsigeps = 0.5*(spsoundi + spsoundj)
       !vsigeps = 0.25*(dustfraci + dustfracj)*(spsoundi + spsoundj)
       !diffeps = alphaB*rhoav1*vsigeps*(dustfraci - dustfracj)*grkern
       !ddustevoldt(i) = ddustevoldt(i) + pmassj*diffeps
       !ddustevoldt(j) = ddustevoldt(j) - pmassi*diffeps
    endif

  end subroutine artificial_dissipation_dust_diffusion

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
!       del2v(i) = del2v(i) + 0.5*pmassj*rhoav1*dx(1)*dr(1)*(-2.*grkern)
!       del2v(j) = del2v(j) + 0.5*pmassi*rhoav1*dx(1)*dr(1)*(-2.*grkern)
       del2v(i) = del2v(i) - pmassj*rhoav1*dvel(1)*grgrw
       del2v(j) = del2v(j) + pmassi*rhoav1*dvel(1)*grgrw
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
    real :: dgradwdhi,dgradwdhj,BdotBextj,term
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
       if (iambipolar > 0) then
       !
       !--curl B has already been calculated, here we use it to compute
       !  the ambipolar diffusion term
       !
          call cross_product3D(curlB(:,i),Bi,vec) ! J x B
          call cross_product3D(vec,Bi,JxBxBi)     ! J x B x B
          call cross_product3D(curlB(:,j),Bj,vec) ! J x B
          call cross_product3D(vec,Bj,JxBxBj)     ! J x B x B
          call cross_product3D(JxBxBi,dr,curlBi)  ! misuse of curlBi variable to store (JxB) x B x \nabla
          call cross_product3D(JxBxBj,dr,curlBj)  ! misuse of curlBi variable to store (JxB) x B x \nabla
          termi = 1./(rhoi*rho_ion*gamma_ambipolar)
          termj = 1./(rhoj*rho_ion*gamma_ambipolar)
          dBdtambi(:) = termi*curlBi(:)*rho21i*grkerni + termj*curlBj(:)*rho21j*grkernj
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

    !--------------------------------------------
    !  real resistivity in induction equation
    !--------------------------------------------
    if (iresist.gt.0) then
       if (imhd.gt.0) then
          dBdtvisc(:) = -2.*etamhd*dB(:)/rij
       else !--vector potential resistivity
          dBdtvisc(:) = 2.*etamhd*dBevol(:)/rij
       endif
       !
       !--add to dB/dt (converted to d(B/rho)/dt later if required)
       !
       dBevoldti(:) = dBevoldti(:) - pmassj*rho1j*dBdtvisc(:)*grkerni
       dBevoldt(:,j) = dBevoldt(:,j) + pmassi*rho1i*dBdtvisc(:)*grkernj
       !
       !  add to thermal energy equation
       !
       if (iener.eq.3 .or. iener.eq.1) then
          stop 'unimplemented iener for physical resistivity'
       elseif (iener.gt.0) then
          term = -etamhd*rho1i*rho1j*dot_product(dB,dB)*grkern/rij
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
    real :: termi,termj,term,dterm,du,si,sj
    real, dimension(ndimV) :: prdustterm

    !
    !--time derivative of dust to gas ratio
    !  (symmetric derivative to conserve dust/gas mass)
    !
    if (use_sqrtdustfrac) then
       si = sqrt(rhoi*dustfraci)
       sj = sqrt(rhoj*dustfracj)
       termi = (1. - dustfraci)*rho1i*projdeltavi*grkerni
       termj = (1. - dustfracj)*rho1j*projdeltavj*grkernj
       term = termi + termj
       ddustevoldt(i) = ddustevoldt(i) - pmassj*sj*term
       ddustevoldt(j) = ddustevoldt(j) + pmassi*si*term
    else
       termi = rhogrhodonrhoi*projdeltavi*rho21i*grkerni
       termj = rhogrhodonrhoj*projdeltavj*rho21j*grkernj
       term  = termi + termj
       ddustevoldt(i) = ddustevoldt(i) - pmassj*term
       ddustevoldt(j) = ddustevoldt(j) + pmassi*term
    endif

    !
    !--time derivative of deltav
    !  (here only bits that involve sums over particles, i.e. not decay term)
    !
    !--high mach number term
    termi    = (rhogasi - rhodusti)*rho1i*deltav2i
    termj    = (rhogasj - rhodustj)*rho1j*deltav2j
    dterm    = 0.5*(termi - termj)
    !
    !--note: need to add term to d/dt(deltav) using the gas-only force 
    !  i.e. before we add the anisotropic pressure term to the forces.
    !  Using forcei and forcej directly means that we automatically include viscosity, MHD etc.
    !
    ddeltavdt(:,i) = ddeltavdt(:,i) + rho1i*pmassj*(dvel(:)*projdeltavi + dterm*dr(:))*grkerni ! add forcei term once finished
    ddeltavdt(:,j) = ddeltavdt(:,j) + rho1j*pmassi*(dvel(:)*projdeltavj + dterm*dr(:))*grkernj - rhoj/rhogasj*forcej(:)

    !
    !--anisotropic pressure term
    !
    prdustterm(:) = rhogrhodonrhoi*deltavi(:)*projdeltavi*rho21i*grkerni &
                  + rhogrhodonrhoj*deltavj(:)*projdeltavj*rho21j*grkernj

    fextrai(:) = fextrai(:) - pmassj*(prdustterm(:))
    fextraj(:) = fextraj(:) + pmassi*(prdustterm(:))

    !
    !--thermal energy equation: add Pg/rhog*div(vgas) and deltav.grad(u) term
    !
    if (iener.gt.0) then
       du = uu(i) - uu(j)
       dudt(i) = dudt(i) + pmassj*(pri*rho1i/rhogasi*projdvgas - rhodusti*rho21i*du*projdeltavi)*grkerni
       dudt(j) = dudt(j) + pmassi*(prj*rho1j/rhogasj*projdvgas - rhodustj*rho21j*du*projdeltavj)*grkernj
    endif

  end subroutine dust_derivs

!----------------------------------------------------------------
!  Derivative of dust fraction and thermal energy
!  for one fluid dust in the terminal velocity approximation
!----------------------------------------------------------------
  subroutine dust_derivs_diffusion
    real :: diffterm, Di, Dj, du, Dav
    real :: tstopi, tstopj, pdvtermi, pdvtermj
    real :: si, sj, rhoav1d

    tstopi = get_tstop(idrag_nature,rhogasi,rhodusti,spsoundi,Kdrag)
    tstopj = get_tstop(idrag_nature,rhogasj,rhodustj,spsoundj,Kdrag)
    
    if (use_sqrtdustfrac) then
       Di = tstopi*rho1i
       Dj = tstopj*rho1j
       Dav = 0.5*(Di + Dj)
       si = sqrt(dustfraci*rhoi)
       sj = sqrt(dustfracj*rhoj)
       !rhoav1d = 0.5*(rho1i + rho1j)
       rhoav1d = 2./(rhoi + rhoj)
       diffterm = rhoav1d*2.*Dav*(pri - prj)*grkern/rij
       ddustevoldt(i) = ddustevoldt(i) - pmassj*sj*diffterm
       ddustevoldt(j) = ddustevoldt(j) + pmassi*si*diffterm
    else
       Di = dustfraci*tstopi
       Dj = dustfracj*tstopj
       !if (Di + Dj > 0.) then
          Dav = 0.5*(Di + Dj)
          !Dav = 2.*Di*Dj/(Di + Dj)
          !Dav = sqrt(0.5*(Di**2 + Dj**2)) !0.5*(Di + Dj) !2.*Di*Dj/(Di + Dj)
       !else
       !   Dav = 0.
       !endif
       diffterm = rho1i*rho1j*2.*Dav*(pri - prj)*grkern/rij
       ddustevoldt(i) = ddustevoldt(i) - pmassj*diffterm
       ddustevoldt(j) = ddustevoldt(j) + pmassi*diffterm
    endif

    !
    !--thermal energy equation: add Pg/rhog*div(vgas) and deltav.grad(u) term
    !
    if (iener.gt.0) then
       du = uu(i) - uu(j)
       pdvtermi = pri*rho1i/rhogasi*pmassj*dvdotr*grkerni
       pdvtermj = prj*rho1j/rhogasj*pmassi*dvdotr*grkernj
       dudt(i) = dudt(i) + pdvtermi - 0.5*pmassj*du*diffterm/(1. - dustfraci)
       dudt(j) = dudt(j) + pdvtermj - 0.5*pmassi*du*diffterm/(1. - dustfracj)
    endif

  end subroutine dust_derivs_diffusion

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
