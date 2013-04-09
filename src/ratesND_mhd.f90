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
 real :: vsig2i,vsig2j,vsigproji,vsigprojj,vsignonlin,vsigu
!! real :: vsigii,vsigjj
 real :: prneti,prnetj,pequil
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
 real :: dusttogasi,dusttogasj
 real :: rhogrhodonrhoi, rhogrhodonrhoj
 real :: deltav2i,deltav2j,rhoonrhog2i
 real :: dtstop,projdeltavi,projdeltavj
 real :: rhogasi,rhodusti,rhogasj,rhodustj,projdvgas
 real, dimension(ndimV) :: vgasi,vgasj,dvgas,fextrai,fextraj
 real :: sum

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
 real :: vsigdtc,zero,fhmax, fonh, forcemag
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
  if (imhd.gt.0) curlB(:,i) = 0.0
  if (allocated(curlBsym)) curlBsym(:,i) = 0.
  if (allocated(divBsym)) divBsym(i) = 0.
  if (allocated(dveldx)) then
     dveldx(:,:,i) = 0.
     dxdx(:,i) = 0.
  endif
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
  if (idust.eq.1) then
     ddusttogasdt(i) = 0.
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
       phii = phi(i)
       phii1 = 1./phii
       sqrtgi = sqrtg(i)
       ! one fluid dust definitions
       if (idust.eq.1) then
          dusttogasi  = dusttogas(i)
          rhogasi        = rhoi/(1. + dusttogasi)
          rhodusti       = rhogasi*dusttogasi
          deltavi(:)     = deltav(:,i)
          deltav2i       = dot_product(deltavi,deltavi)
          rhogrhodonrhoi = rhoi*dusttogasi/(1. + dusttogasi)**2
          vgasi(:)       = veli(:) - dusttogasi/(1. + dusttogasi)*deltavi(:)
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
          alphaBi = alpha(3,i)
          if (imhd.lt.0) Bevoli(:) = Bevol(:,i)
          if (iresist.eq.1) then
             etai = etamhd
          elseif (iresist.eq.3) then
             etai = etafunc(x(1,i),etamhd)
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
              if (rij.le.1.e-8 .and. itype(j).eq.itypei) then
                 nclumped = nclumped + 1
                 if (rij.lt.tiny(rij)) then
                    write(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                    call quit
                 endif
              endif
              dr(1:ndim) = dx(1:ndim)/rij      ! unit vector           
              !do idim=1,ndim
              !   dri(idim) = dot_product(1./gradmatrix(idim,1:ndim,i),dr(1:ndim))
              !   drj(idim) = dot_product(1./gradmatrix(idim,1:ndim,j),dr(1:ndim))
              !enddo
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
       force(:,i) = force(:,i) + forcei(:) + fextrai(:)
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
 do i=1,npart
 
    rhoi  = rho(i)
    rho1i = 1./rhoi
    if (ivisc.gt.0) then
       sum = sum + pmass(i)*dot_product(vel(:,i),force(:,i))
    endif
!
!--Dust
!
    if (idust.eq.2 .and. idrag_nature.ne.0 .and. Kdrag.gt.0. .and. ismooth.lt.1) then
       !
       !--two fluid dust: calculate drag timestep
       !
       dtdrag = min(dtdrag,rhoi/Kdrag)
    elseif (idust.eq.1) then
       !print *,'dtdrag = ',min(dtdrag,rhoi/Kdrag)
       !dtdrag = min(dtdrag,0.25*rhoi/Kdrag)
       !------------------
       !  one fluid dust
       !------------------
       !
       !--d/dt(rhod/rhog): multiply by terms out the front of the SPH sum
       !
       rhoonrhog2i = (1. + dusttogas(i))**2  ! (rho/rhog)**2
       ddusttogasdt(i) = rhoonrhog2i*ddusttogasdt(i)

       dusttogasi = dusttogas(i)
       rhogasi = rhoi/(1. + dusttogasi)
       rhodusti = rhogasi*dusttogasi
       
       !
       !--d/dt(deltav)  : add terms that do not involve sums over particles
       !
       if (dusttogas(i).gt.0.) then
          dtstop   = Kdrag*(1. + dusttogas(i))**2/(rhoi*dusttogas(i))  ! 1/tstop = K*rho/(rhod*rhog)
          !ddeltavdt(:,i) = ddeltavdt(:,i) - deltav(:,i)*dtstop
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
       !
       !--DEBUGGING: CHECK ENERGY CONSERVATION
       !  BY ADDING TERMS (SHOULD GIVE ZERO)
       !
       sum = sum + pmass(i)*(dot_product(vel(:,i),force(:,i)) &
                 + rhogasi*rhodusti*rho1i**2*dot_product(deltav(:,i),ddeltavdt(:,i)) &
                 + rhogasi**2*rho1i**2*((1. - dusttogasi)/(1. + dusttogasi)* &
                   0.5*dot_product(deltav(:,i),deltav(:,i)) - uu(i))*ddusttogasdt(i) &
                 + rhogasi*rho1i*dudt(i))
    endif

!
!--compute JxB force from the Euler potentials (if using second derivs)
!  also divide by rho for J and div B calculated the "normal" way
!
    if (imhd.ne.0) then
       if (imagforce.eq.4) then
          call cross_product3D(curlB(:,i),Bfield(:,i),fmagi(:)) ! J x B
          force(:,i) = force(:,i) + fmagi(:)*rho1i  ! (J x B)/rho
          fmag(:,i) = fmagi(:)*rho1i
       elseif (imhd.gt.0) then ! not for vector potential
          curlB(:,i) = curlB(:,i)*rho1i
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
                            ndim,ndimV,vel(1:ndimV,i),hh(i),spsound(i),itype(i), &
                            rhoi)
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
!--calculate resistive timestep (bootstrap onto force timestep)
!
    if (iresist.gt.0 .and. iresist.ne.2 .and. etamhd.gt.tiny(etamhd)) then
       fhmax = max(fhmax,etai/(hh(i)**2))
    endif

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
    elseif (iener.gt.0 .and. iav.ge.0 .and. iener.ne.10 .and. iener.ne.11 .and. idust.ne.1) then
       if (damp.lt.tiny(0.)) dudt(i) = dudt(i) + pr(i)*rho1i**2*drhodt(i)    
       dendt(i) = dudt(i)
    else
       dendt(i) = dudt(i)
    endif
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
 enddo
 if (idust.eq.1 .and. abs(sum).gt.1.e-9) print*,' SUM (should be zero if conserving energy) = ',sum
 if (ivisc.gt.0) print*,' dEk/dt = ',sum
 
 if (sqrt(dot_product(fmean,fmean)).gt.1.e-8 .and. mod(nsteps,100).eq.0) print*,'WARNING: fmean = ',fmean(:)
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
  subroutine drag_forces
!--DIRTY declarations for the hack  
    use kernels, only:interpolate_kerneldrag,interpolate_kernel  
    use options, only:idrag_nature,idrag_structure,Kdrag   
    implicit none
    integer :: itypej
    logical :: iskip_drag
    real    :: coeff_gei_1,coeff_dq_1,coeff_dq_4
    real    :: dv2,vij,V0,f,dragcoeff,dragterm,dragterm_en,dragcheck
    real    :: s2_over_m,spsoundgas
    real    :: wabj,wab,hfacwabj,rhoiplusj,rhoiplusj2
    real    :: dt1,dt12
    !real    :: gkeri,gkerj
    real, dimension(ndimV) :: drdrag
    real, parameter :: pow_drag_exp = 0.4
    real, parameter :: a2_set       = 0.5
    real, parameter :: a3_set       = 0.5
    real, parameter :: pi  = 3.141592653589
!
!--setup the parameters
!
    iskip_drag  = .false.
    coeff_gei_1 = 4./3.*sqrt(8.*pi/gamma)
    coeff_dq_1  = sqrt(0.5*gamma)
    coeff_dq_4  = 9.*pi*gamma/128.
    s2_over_m   = 1./coeff_gei_1
    f           = 0.
    V0          = 0.
    dragcoeff   = 0.
    dragcheck   = 0.
!    dt1         = 1./tout
          dt1         = 1./(1.8716d-2)
    dt12        = dt1*dt1
!
!--get informations on the differential velocities
!
    velj(:) = vel(:,j)
    dvel(:) = veli(:) - velj(:)
    dv2     = dot_product(dvel,dvel)
    
    if (rij.le.epsilon(rij)) then !two particles at the same position
       if (dv2.le.epsilon(dv2)) then !no differential velocity => no drag
           drdrag(1:ndim) = 0.
           iskip_drag      = .true.
       else ! Change dr so that the drag is colinear to the differential velocity
           vij            = sqrt(dv2)
           drdrag(:) = dvel(:)/vij
       endif
    else ! dr = drdrag
      drdrag(:)       = dr(:)
    endif

!---start the drag calculation    
    !if (.not.iskip_drag) then
!
!--get the j particle extra properties
!
!--Hack special SI    
    itypej     = itype(j)
    pmassj     = pmass(j)
    rhoj       = rho(j)
    rhoij      = rhoi*rhoj
    rhoiplusj  = rhoi+rhoj
    rhoiplusj2 = rhoiplusj*rhoiplusj

    if (itypei.eq.itypegas) then
       spsoundgas = spsound(i)
    else
       spsoundgas = spsound(j)
    endif
!
!--calculate the kernel(s)
!
    hfacwabj = (1./hh(j)**ndim)
!-- HACK
!    call interpolate_kernel(q2i,wabi,gkeri)
!--OK
    call interpolate_kerneldrag(q2i,wabi)
    wabi     = wabi*hfacwabi
!--DIRTY HACK
!   call interpolate_kernel(q2j,wabj,gkerj) 
!-OK
   call interpolate_kerneldrag(q2j,wabj)
    wabj     = wabj*hfacwabj
    !wab = 0.5*(wabi + wabj)
    if (itypei.eq.itypegas) then
       wab = wabi
    else
       wab = wabj
    endif   
!
!--calculate the quantities needed for the drag force
!
    V0 = dot_product(dvel,drdrag)
    
    select case(idrag_nature)
        case(1) !--constant drag
           dragcoeff = Kdrag/rhoij
        case(2) !--Epstein regime
           select case(idrag_structure)
              case(1,4,5) !--linear, third order, PM expression
                 dragcoeff = coeff_gei_1*s2_over_m*spsoundgas
              case(3)
                 dragcoeff = pi*s2_over_m
              case default
                 print*,'ERROR drag calculation: wrong idrag_structure'
           end select
    end select
    
    select case(idrag_structure)
       case(1) !--linear regime
          f = 1.
       case(2) !--power law
          f = dv2**(0.5*pow_drag_exp)
       case(3) !--quadratic
          f = sqrt(dv2)
       case(4) !--cubic expansion
          select case(idrag_nature)
             case(1) !--forced drag
                f = 1. + a3_set*dv2
             case(2) !--Epstein
                f = 1. + 0.2*coeff_dq_1*coeff_dq_1*dv2/(spsoundgas*spsoundgas)
             case default
                print*,'ERROR drag calculation cubic: wrong drag nature' 
          end select
       case(5) !--PM expression
          select case(idrag_nature)
             case(1) !--forced drag
                f = sqrt(1. + a2_set*dv2)
             case(2) !--Epstein
                f = sqrt(1. + coeff_dq_4*dv2/(spsoundgas*spsoundgas))
             case default
                print*,'ERROR drag calculation PM: wrong drag nature'          
          end select
       case default
          print*,'this value for idrag_structure does not exist'         
    end select
!
!--update the force and the energy
!
 
!--DIRTY HACK
!    dragterm    = wab*dragcoeff*f
!--OK
    dragterm    = ndim*wab*dragcoeff*f*V0        
    dragterm_en = dragterm*V0
!--DIRTY HACK    
!    forcei(:)   = forcei(:)  - dragterm*pmassj*dvel(:)
!    force(:,j)  = force(:,j) + dragterm*pmassi*dvel(:)    
!-OK
    forcei(:)   = forcei(:)  - dragterm*pmassj*drdrag(:)   
!-OK
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
    use cons2prim, only:specialrelativity
    implicit none
    real :: prstar,prstari,prstarj,vstar
    real :: dvsigdtc
    real :: hav,hav1,h21,q2
    real :: hfacwab,hfacwabj,hfacgrkern,hfacgrkernj
    real :: wabalti,wabaltj,wabalt
    real :: altrhoi,altrhoj,gammastar,vperp2
    real :: enthalpi,enthalpj,term,denom,term1

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
    if (idust.eq.1) then
       dusttogasj  = dusttogas(j)
       rhogasj        = rhoj/(1. + dusttogasj)
       rhodustj       = rhogasj*dusttogasj
       deltavj(:)     = deltav(:,j)
       deltav2j       = dot_product(deltavj,deltavj)
       rhogrhodonrhoj = rhoj*dusttogasj/(1. + dusttogasj)**2
       vgasj(:)       = velj(:) - dusttogasj/(1. + dusttogasj)*deltavj(:)
       dvgas(:)       = vgasi(:) - vgasj(:)
       projdvgas      = dot_product(dvgas,dr)
    else
       projdvgas      = dvdotr
       deltav2j       = 0.
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
       if (iresist.eq.1) then
          etaj = etamhd
       elseif (iresist.eq.3) then
          etaj = etafunc(x(1,j),etamhd)
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
          vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
          vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))
          vsignonlin = 0.5*(vsigi + vsigj) + abs(dvdotr)
          !vsignonlin = sqrt(dot_product(dvel - dvdotr,dvel - dvdotr))
          !vsignonlin = 0.5*(sqrt(valfven2i) + sqrt(valfven2j))
       else
          vsigi = spsoundi
          vsigj = spsoundj
          vsignonlin = 0.
       endif

       vsig = 0.5*(max(vsigi + vsigj - beta*dvdotr,0.0)) ! also used where dvdotr>0 in MHD
       !vsig = 0.5*(vsigi + vsigj + beta*abs(dvdotr))
       !vsignonlin = 0.5*max(vsigi + vsigj - 4.0*dvdotr,0.0) !!!*(1./(1.+exp(1000.*dvdotr/vsig)))
       !vsignonlin = max(-dvdotr,0.0) !!!*(1./(1.+exp(1000.*dvdotr/vsig)))

       vsigu = sqrt(abs(prneti-prnetj)*rhoav1)
       !vsigu = sqrt(0.5*(psi(i) + psi(j)))
       !vsigu = abs(dvdotr)

       ! vsigdtc is the signal velocity used in the timestep control
       vsigdtc = max(0.5*(vsigi + vsigj + beta*abs(dvdotr)),vsignonlin)
       if (idust.eq.1) then
          vsigdtc = vsigdtc + sqrt(deltav2i + deltav2j)
       endif
    
    endif

    dvsigdtc = 1./vsigdtc
    vsigmax = max(vsigmax,vsigdtc)
    !
    !--time step control (courant and viscous)
    !
    if (vsigdtc.gt.zero) dtcourant = min(dtcourant,min(hi*dvsigdtc,hj*dvsigdtc))
    
!----------------------------------------------------------------------------
!  artificial dissipation terms
!----------------------------------------------------------------------------

    if (iav.gt.0) then
       if (specialrelativity) then
          call artificial_dissipation_sr
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
          !!prterm = (pri - prj)/rhoij*grkern

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

!    forcei(:) = forcei(:) + pmassj*(pr(i)/rho(i)*dri(:)*grkerni &
!                            + pr(j)/rho(j)*drj(:)*grkernj)
!    forcej(:) = forcej(:) - pmassi*(pr(i)/rho(i)*dri(:)*grkerni &
!                            + pr(j)/rho(j)*drj(:)*grkernj)

!------------------------------------------------------------------------
!  Self-gravity term
!------------------------------------------------------------------------
     
!------------------------------------------------------------------------
!  Lorentz force and time derivative of B terms
!------------------------------------------------------------------------

    if (imhd.ne.0) call mhd_terms

!------------------------------------------------------------------------
!  Time derivative of dust to gas ratio and deltav for one fluid dust
!------------------------------------------------------------------------

    if (idust.eq.1) call dust_derivs
!
!   Add contributions to j from mhd terms and dissipation here
!   to avoid repeated memory access
!
!   here forcej is the gas-only forces, fextraj are forces that act on the total fluid
    force(:,j) = force(:,j) + forcej(:) + fextraj(:)

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
    real :: vissv,vissB,vissu,vissdust
    real :: term,dpmomdotr
    real :: termnonlin,termu,termv,termi,termj
    real :: rhogonrhoi,rhogonrhoj,rhogonrhoav,rhoav1g
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
    elseif (idust.eq.1) then    
       dpmomdotr = -projdvgas
    else
       dpmomdotr = -dvdotr
    endif
    !--used for viscosity
    term = vsig*rhoav1*grkern
    termv = term
    
    if (idust.eq.1) then
    !
    !--for one fluid dust, viscosity term must be multiplied by
    !  rhog/rho, but done in a symmetric way to conserve momentum
    !  rhog/rho = 1/(1 + dusttogas)
    !
       rhogonrhoi = 1./(1. + dusttogasi)
       rhogonrhoj = 1./(1. + dusttogasj)
       rhogonrhoav = 0.5*(rhogonrhoi + rhogonrhoj)
       termv = termv*rhogonrhoav
       if (iav.ge.3) stop 'av terms not implemented for iav>=3 and one-fluid dust'
       if (iener.eq.3) stop 'av terms not implemented for iener=3 and one-fluid dust'
    endif

    !--used for thermal conductivity
    termu = vsigu*rhoav1*grkern
    !termu = vsig*rho1i*rho1j*grkern
    !termu = vsig*rhoav1*grkern

    !--used for resistivity
    termnonlin = vsignonlin*rhoav1*grkern
!    termnonlin = vsigi*alphaBi*vsigj*alpha(3,j)/ &
!                (vsigi*alphaBi + vsigj*alpha(3,j))*rhoav1*grkern
!    termnonlin = vsigi*alphaBi*vsigj*alpha(3,j)*grkerni*grkernj/ &
!                (rhoi*vsigi*alphaBi*grkerni + rhoj*vsigj*alpha(3,j)*grkernj)
    !termnonlin = vsignonlin*rhoav1*grkern

    
    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------
    
    if (dvdotr.lt.0 .and. iav.le.2) then            
       visc = alphaav*termv*dpmomdotr     ! viss=abs(dvdotr) defined in rates
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
    elseif (iav.ge.3) then ! using total energy, for approaching and receding
       visc = alphaav*termv
       !print*,'visc = ',i,j,vsig*alphaav*hh(i)
       visc = visc*0.5*(hh(i) + hh(j))/rij  ! use this line to multiply viscosity by h/r
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
          dBdtvisc(:) = alphaB*termnonlin*Bvisc(:)

          !
          !--add to d(B/rho)/dt (converted to dB/dt later if required)
          !
          dBevoldti(:) = dBevoldti(:) + rhoi*pmassj*dBdtvisc(:)               
          dBevoldt(:,j) = dBevoldt(:,j) - rhoj*pmassi*dBdtvisc(:)

       elseif (iav.eq.1) then !--wrong vector potential resistivity
          dBdtvisc(:) = alphaB*termnonlin*dBevol(:)
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
       if (dvdotr.lt.0 .and. iav.le.2) then
          v2i = dot_product(veli,dr)**2      ! energy along line
          v2j = dot_product(velj,dr)**2      ! of sight
          qdiff = qdiff + alphaav*0.5*(v2i-v2j)
       elseif (iav.ge.3) then
          v2i = dot_product(veli,veli)      ! total energy
          v2j = dot_product(velj,velj)
          qdiff = qdiff + alphaav*0.5*(v2i-v2j)
       endif
       !
       !  thermal energy terms - applied everywhere
       !
       qdiff = qdiff + alphau*(uu(i)-uu(j))
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
          qdiff = qdiff + alphaB*0.5*(B2i-B2j)*rhoav1
       elseif (imhd.lt.0) then
          stop 'mhd dissipation not implemented with total energy equation for vector potential'
       endif
       !
       !  add to total energy equation
       !
       dendt(i) = dendt(i) + pmassj*term*qdiff
       dendt(j) = dendt(j) - pmassi*term*qdiff

    !--------------------------------------------------
    !   thermal energy equation
    !--------------------------------------------------       
    elseif (iener.gt.0) then
       !
       !  kinetic energy terms
       !
       if (dvdotr.lt.0 .and. iav.le.2) then
          if (idust.eq.1) then
             vissv = -alphaav*0.5*projdvgas**2          
          else
             vissv = -alphaav*0.5*(dot_product(veli,dr) - dot_product(velj,dr))**2
          endif
       elseif (iav.ge.3) then
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
          if (idust.eq.1) then
             dudt(i) = dudt(i) + pmassj*(term*(vissv) + termu*rhoi/rhogasi*vissu + termnonlin*(vissB))
             dudt(j) = dudt(j) + pmassi*(term*(vissv) - termu*rhoj/rhogasj*vissu + termnonlin*(vissB))          
          else
             dudt(i) = dudt(i) + pmassj*(term*(vissv) + termu*vissu + termnonlin*(vissB))
             dudt(j) = dudt(j) + pmassi*(term*(vissv) - termu*vissu + termnonlin*(vissB))
          endif

          !--entropy dissipation
          if (iener.eq.1) then
             if (alphau.gt.0.) stop 'thermal conduction in entropy equation is wrong'
             vissu = alphau*(en(i)-en(j))
             dendt(i) = dendt(i) + pmassj*(termu*(vissu))
             dendt(j) = dendt(j) + pmassi*(termu*(-vissu))
          endif
       endif

       !
       !--dissipation term in dust-to-gas ratio
       !
       if (idust.eq.1) then
!    termu = vsigu*rhoav1*grkern
          rhoav1g = 1./(0.5*(rhogasi + rhogasj))
          vissdust = vsig*rhoav1*grkern*alphaB*(rhodusti - rhodustj)*rhoav1g
          ddusttogasdt(i) = ddusttogasdt(i) + pmassj*(vissdust)
          ddusttogasdt(j) = ddusttogasdt(j) - pmassi*(vissdust)
          dudt(i) = dudt(i) + uu(i)*rhoi/rhogasi*pmassj*(vissdust)
          dudt(j) = dudt(j) - uu(j)*rhoj/rhogasj*pmassi*(vissdust)
          termi = 0.5*rhoi*rhoi/(rhogasi*rhodusti)*(1. - dusttogasi)/(1. + dusttogasi)
          termj = 0.5*rhoj*rhoj/(rhogasj*rhodustj)*(1. - dusttogasj)/(1. + dusttogasj)
          ddeltavdt(:,i) = ddeltavdt(:,i) - deltavi(:)*termi*pmassj*(vissdust)
          ddeltavdt(:,j) = ddeltavdt(:,j) + deltavj(:)*termj*pmassi*(vissdust)
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
    real :: vissv,term,visc
    real :: alphaav,alphau,alphaB
    !
    !--definitions
    !      
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    alphaB = 0.5*(alphaBi + alpha(3,j))

    vissv = -alphaav*0.5*(projdvgas)**2
    term = vsig*rhoav1*grkern

    if (dvdotr.lt.0 .and. iav.le.2) then            
       visc = alphaav*term*abs(projdvgas)
       forcei(:) = forcei(:) - pmassj*visc*dr(:)
       forcej(:) = forcej(:) + pmassi*visc*dr(:)
    endif

  end subroutine artificial_dissipation_phantom

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
    real :: termnonlin,lorentzfactorstari,lorentzfactorstarj,estari,estarj
    real :: lorentzi,lorentzj
    real,dimension(ndimV) :: pmomstari,pmomstarj
    !
    !--definitions
    !       
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    vsigav = max(alphaav,alphau)*vsig
    !!rhoav1 = 2./(rhoi + rhoj)
!    if (geom(1:4).ne.'cart') then
!    dpmomdotr = abs(dot_product(pmom(:,i)-pmom(:,j),dr(:)))
!    else
!       dpmomdotr = -dvdotr
!    endif
    term = vsig*rhoav1*grkern
    termnonlin = vsignonlin*rhoav1*grkern

    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------
    
    lorentzfactorstari=1.0/(sqrt(1.0-(projvi**2)))
    lorentzfactorstarj=1.0/(sqrt(1.0-(projvj**2)))
    
    dens1j = 1./dens(j)
    estari=lorentzfactorstari*(1.0+alphau*(uu(i))+0.*(pri*dens1i))-0.*(pri*rho1i)
    estarj=lorentzfactorstarj*(1.0+alphau*(uu(j))+0.*(prj*dens1j))-0.*(prj*rho1j)
    
!    pmomstari=veli*(estari+0.*(pri*rho1i))
!    pmomstarj=velj*(estarj+0.*(prj*rho1j))
    
    pmomstari = veli*lorentzfactorstari
    pmomstarj = velj*lorentzfactorstarj

    dpmomdotr=abs(dot_product(pmomstari-pmomstarj,dr))
    
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
          qdiff = qdiff + alphaav*(estari-estarj)
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
          qdiffi=(estari-estarj)+(projvi*dpmomdotr)
          qdiffj=(estarj-estari)-(projvj*dpmomdotr)
       endif
       
       if (v2i.le.1.0) lorentzi=1.0/(sqrt(1.0-(v2i)))
       if (v2j.le.1.0) lorentzj=1.0/(sqrt(1.0-(v2j)))
       !
       !  add to entropy equation
       !
       !if (qdiffi.gt.1.e-8 .or. qdiffj.gt.1.e-8) then
       !   print*,'error, negative entropy generation ',qdiffi,qdiffj
          !stop
       !endif
       dudt(i) = (dudt(i) + alphaav*lorentzi*pmassj*term*qdiffi)
       dudt(j) = (dudt(j) + alphaav*lorentzj*pmassi*term*qdiffj)

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
    real, dimension(ndimV) :: dBdtvisc,curlBj,curlBterm
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
    if (imagforce.eq.4) then
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
       !
       !--calculate curl B for current (only if field is 3D)
       !  this is used in the switch for the artificial resistivity term
       !
       if (imhd.gt.0) then ! not for vector potential
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
          - phii_on_phij*pmassj*(dvel(:)*projBrhoi)*grkerni
       dBevoldt(:,j) = dBevoldt(:,j)          &
          - phij_on_phii*pmassi*(dvel(:)*projBrhoj)*grkernj
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
    if (iresist.gt.0 .and. iresist.ne.2) then
       !etaij = 0.5*(etai + etaj)
       etaij = 2.*etai*etaj/(etai + etaj)
       if (imhd.gt.0) then
          dBdtvisc(:) = -2.*etaij*dB(:)/rij
       else !--vector potential resistivity
          dBdtvisc(:) = 2.*etaij*dBevol(:)/rij
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
!  for one fluid dust
!----------------------------------------------------------------
  subroutine dust_derivs
    implicit none
    real :: termi,termj,term,dterm,du
    real, dimension(ndimV) :: prdustterm

    !
    !--time derivative of dust to gas ratio
    !  (symmetric derivative to conserve dust/gas mass)
    !
    projdeltavi = dot_product(deltavi,dr)
    projdeltavj = dot_product(deltavj,dr)
    termi = rhogrhodonrhoi*projdeltavi*rho21i*grkerni
    termj = rhogrhodonrhoj*projdeltavj*rho21j*grkernj
    term  = termi + termj
    
    ddusttogasdt(i) = ddusttogasdt(i) - pmassj*term
    ddusttogasdt(j) = ddusttogasdt(j) + pmassi*term

    !
    !--time derivative of deltav
    !  (here only bits that involve sums over particles, i.e. not decay term)
    !
    !--high mach number term
    rhogasj  = rhoj/(1. + dusttogasj)
    rhodustj = rhogasj*dusttogasj
    termi    = (1. - dusttogasi)/(1. + dusttogasi)*deltav2i
    termj    = (1. - dusttogasj)/(1. + dusttogasj)*deltav2j
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
