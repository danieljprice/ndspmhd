!!--------------------------------------------------------------------
!! Computes the rates of change of the conserved variables
!! (forces, energy etc)
!! This is the core of the SPH algorithm
!!--------------------------------------------------------------------

subroutine get_rates
! USE dimen_mhd
 use debug, only:trace
 use loguns, only:iprint
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
 use grutils, only:metric_diag,dot_product_gr
 !!use matrixcorr
!
!--define local variables
!
 implicit none
 integer :: i,j,n
 integer :: icell,iprev,nneigh,nlistdim
 integer, allocatable, dimension(:) :: listneigh
 integer :: idone
!
!  (particle properties - local copies and composites)
!
 real :: rij,rij2
 real :: rhoi,rho1i,rho2i,rho21i,rhoj,rho1j,rho2j,rho21j,rhoav1,rhoij
 real :: pmassi,pmassj
 real :: Prho2i,Prho2j,prterm,pri,prj
 real :: hi,hi1,hj,hj1,hi21,hj21
 real :: hfacwabi,hfacgrkerni
 real, dimension(ndim) :: dx, fexternal  !!,dri,drj
!
!--gr terms
!
 real :: sqrtgi,sqrtgj
 real, dimension(ndimV) :: gdiagi,gdiagj
!
!  (velocity)
!      
 real, dimension(ndimV) :: veli,velj,dvel
 real, dimension(ndimV) :: dr
 real :: dvdotr
!
!  (mhd)
!     
 real, dimension(ndimB) :: Brhoi,Brhoj,Bi,Bj,dB
 real, dimension(ndimB) :: faniso,fmagi,fmagj
 real, dimension(ndimB) :: curlBi
 real :: fiso,B2i,B2j
 real :: valfven2i,valfven2j
 real :: BidotdB,BjdotdB,Brho2i,Brho2j
 real :: projBrhoi,projBrhoj,projBi,projBj,projdB,projBconst
!
!  (artificial viscosity quantities)
!      
 real :: vsig,vsigi,vsigj
 real :: spsoundi,spsoundj,alphai,alphaui,alphaBi
!! real :: rhoi5,rhoj5
 real :: vsig2i,vsig2j,vsigproji,vsigprojj !!,vsignonlin
!! real :: vsigii,vsigjj
!
!  (av switch)
!
 real :: source,tdecay1,sourcedivB,sourceJ,sourceB,sourceu
 real :: graduterm, graddivvmag
 real, dimension(:), allocatable :: del2u
 real, dimension(:,:), allocatable :: graddivv
!
!  (alternative forms)
!
 real, dimension(:), allocatable :: phi
 real :: phii,phii1,phii_on_phij,phij_on_phii
!
!  (kernel related quantities)
!
 real :: q2i,q2j
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj
 real :: gradhi,gradhni
!
!  (time step criteria)
!      
 real :: vsigdtc,zero,fhmax, fonh, forcemag
 integer :: ierr
!
!--div B correction
! 
 real :: gradpsiterm,vsig2,vsigmax !!,dtcourant2
 real :: stressterm, stressmax
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
 listneigh = 0
!
!--initialise quantities
!      
 dtcourant = 1.e6  
 zero = 1.e-10
 vsigmax = 0.
 dr(:) = 0.
 
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
  curlB(:,i) = 0.0
  xsphterm(:,i) = 0.0
  del2u(i) = 0.0
  graddivv(:,i) = 0.0
 enddo

!
!--calculate maximum neg stress for instability correction
!  
 stressmax = 0.
 if (imhd.ne.0 .and. imagforce.eq.2) then
    do i=1,ntotal
       call metric_diag(x(:,i),gdiagi(:),sqrtgi,ndim,ndimV,geom)
       B2i = dot_product_gr(Bfield(:,i),Bfield(:,i),gdiagi(:))
       stressterm = max(0.5*B2i - pr(i),0.)
       stressmax = max(stressterm,stressmax)
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
 if (iprterm.lt.0 .and. iav.eq.0 .and. imhd.eq.0 .and. iener.lt.3) then
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
       rhoi = rho(i)
       rho2i = rhoi*rhoi
!!       rhoi5 = sqrt(rhoi)
       rho1i = 1./rhoi
       rho21i = rho1i*rho1i
       pri = pr(i)
       pmassi = pmass(i)
       Prho2i = pr(i)*rho21i
       spsoundi = spsound(i)
       veli(:) = vel(:,i)
       alphai = alpha(1,i)
       alphaui = alpha(2,i)
       phii = phi(i)
       phii1 = 1./phii
       sqrtgi = sqrtg(i)
       ! mhd definitions
       if (imhd.ne.0) then
          Bi(:) = Bfield(:,i)
          Brhoi(:) = Bi(:)*rho1i
!          if (geom(1:6).ne.'cartes') then
             call metric_diag(x(:,i),gdiagi(:),sqrtgi,ndim,ndimV,geom)
             B2i = dot_product_gr(Bi,Bi,gdiagi)
!          else
!             B2i = dot_product(Bi,Bi)
!          endif
          Brho2i = B2i*rho21i
          valfven2i = B2i*rho1i
          alphaBi = alpha(3,i)
       endif
       gradhi = gradh(i)
       gradhni = gradhn(i)
       hi = hh(i)
       if (hi.le.0.) then
          write(iprint,*) ' rates: h <= 0 particle',i,hi
        call quit
       endif
       hi1 = 1./hi
       hi21 = hi1*hi1
       hfacwabi = hi1**ndim
       hfacgrkerni = hfacwabi*hi1
!
!--for each particle in the current cell, loop over its neighbours
!
       loop_over_neighbours: do n = idone+1,nneigh
       
           j = listneigh(n)
           if ((j.ne.i).and..not.(j.gt.npart .and. i.gt.npart)) then            ! don't count particle with itself
           dx(:) = x(:,i) - x(:,j)
           !print*,' ... neighbour, h=',j,hh(j),rho(j),x(:,j)
           hj = hh(j)
           hj1 = 1./hj
           hj21 = hj1*hj1
     
           rij2 = dot_product(dx,dx)
           q2i = rij2*hi21
           q2j = rij2*hj21

             !----------------------------------------------------------------------------
             !  do pairwise interaction if either particle is within range of the other
             !----------------------------------------------------------------------------
           if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then  ! if < 2h
              rij = sqrt(rij2)
              if (rij.le.1.e-10) then
                 write(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                 call quit
              endif      
              dr(1:ndim) = dx(1:ndim)/rij      ! unit vector
              !do idim=1,ndim
              !   dri(idim) = dot_product(1./gradmatrix(idim,1:ndim,i),dr(1:ndim))
              !   drj(idim) = dot_product(1./gradmatrix(idim,1:ndim,j),dr(1:ndim))
              !enddo
              call rates_core
           else      ! if outside 2h
!              PRINT*,'outside 2h, not calculated, r/h=',sqrt(q2i),sqrt(q2j)
           endif      ! r < 2h or not

         endif            ! j.ne.i
        
        enddo loop_over_neighbours

       iprev = i
       if (iprev.ne.-1) i = ll(i)            ! possibly should be only IF (iprev.NE.-1)
    
    enddo loop_over_cell_particles
            
 enddo loop_over_cells
 
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

 do i=1,npart
 
    rho1i = 1./rho(i)
!
!--compute JxB force from the Euler potentials (if using second derivs)
!  also divide by rho for J and div B calculated the "normal" way
!
    if (imhd.ne.0) then
       if (imagforce.eq.4) then
          call cross_product3D(curlB(:,i),Bfield(:,i),fmagi(:)) ! J x B
          force(:,i) = force(:,i) + fmagi(:)*rho1i  ! (J x B)/rho
          fmag(:,i) = fmagi(:)*rho1i
       else
          curlB(:,i) = curlB(:,i)*rho1i
       endif
       divB(i) = divB(i)*rho1i
    endif
!
!--add external (body) forces
!
    if (iexternal_force.ne.0) then
       call external_forces(iexternal_force,x(1:ndim,i),fexternal(1:ndim),ndim)
       force(1:ndim,i) = force(1:ndim,i) + fexternal(1:ndim)
    endif
!
!--add source terms (derivatives of metric) to momentum equation
!
    if (allocated(sourceterms)) force(:,i) = force(:,i) + sourceterms(:,i)
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
       dendt(i) = (gamma-1.)/dens(i)**(gamma-1.)*dudt(i)      
    elseif (iener.gt.0 .and. iav.ge.0) then
       dudt(i) = dudt(i) + pr(i)*rho1i**2*drhodt(i)    
       dendt(i) = dudt(i)
    else
       dendt(i) = dudt(i)
    endif
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
    elseif (imhd.lt.0) then ! vector potential evolution
       dBevoldt(:,i) = dBevoldt(:,i)*rho1i
    else
       dBevoldt(:,i) = 0.
    endif
!
!--calculate maximum force/h for the timestep condition
!  also check for errors in the force
!
    if ( any(force(:,i).gt.1.e8)) then
       write(iprint,*) 'rates: force ridiculous ',force(:,i),' particle ',i
       call quit
    endif
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
          sourceu = hh(i)*abs(del2u(i))/sqrt(uu(i))
          daldt(2,i) = (alphaumin - alpha(2,i))*tdecay1 + sourceu
       endif
       !
       !--artificial resistivity parameter if iavlim > 10
       !
       if (iavlim(3).ne.0 .and. imhd.ne.0) then
          !--calculate source term for the resistivity parameter
          sourceJ = SQRT(DOT_PRODUCT(curlB(:,i),curlB(:,i))*rho1i)
          sourcedivB = 10.*abs(divB(i))*SQRT(rho1i)
          sourceB = MAX(sourceJ,sourcedivB)
          if (iavlim(3).eq.2) sourceB = sourceB*(2.0-alpha(3,i))
          daldt(3,i) = (alphaBmin - alpha(3,i))*tdecay1 + sourceB
       endif
    endif
!
!--calculate time derivative of divergence correction parameter psi
!      
    select case(idivBzero)
       case(2:7)
          dpsidt(i) = -vsig2max*divB(i) - psidecayfact*psi(i)*vsigmax/hh(i)         
       case DEFAULT
          dpsidt(i) = 0.
    end select
 enddo
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
    if (itype(i).eq.1 .or. i.gt.npart) then
       force(:,i) = 0.0
       drhodt(i) = 0.0
       dudt(i) = 0.0
       dendt(i) = 0.0
       dBevoldt(:,i) = 0.0
       daldt(:,i) = 0.
       dpsidt(i) = 0.0
       fmag(:,i) = 0.0
       divB(i) = 0.0
       curlB(:,i) = 0.0
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
!--------------------------------------------------------------------------------------
! This is the interaction between the particle pairs
! Note that local variables are used from get_rates
!--------------------------------------------------------------------------------------
  subroutine rates_core
    use kernels, only:interpolate_kernel,interpolate_kernels
    implicit none
    real :: prstar, vstar
    real :: projvi, projvj, dvsigdtc
    real :: hav,hav1,h21,q2
    real :: hfacwab,hfacwabj,hfacgrkern,hfacgrkernj
    real :: wabalti,wabaltj,wabalt,grkernalti,grkernaltj

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
      call interpolate_kernels(q2i,wabi,grkerni,wabalti,grkernalti)
      wabi = wabi*hfacwabi
      wabalti = wabalti*hfacwabi
      grkerni = grkerni*hfacgrkerni
      grkernalti = grkernalti*hfacgrkerni
      !  (using hj)
      hfacwabj = hj1**ndim
      hfacgrkernj = hfacwabj*hj1
      call interpolate_kernels(q2j,wabj,grkernj,wabaltj,grkernaltj)
      wabj = wabj*hfacwabj
      wabaltj = wabaltj*hfacwabj
      grkernj = grkernj*hfacgrkernj
      grkernaltj = grkernaltj*hfacgrkernj
      !  (calculate average)              
      wab = 0.5*(wabi + wabj)
      wabalt = 0.5*(wabalti + wabaltj)
      !  (grad h terms)  
      if (ikernav.eq.3) then  ! if using grad h corrections
         !--use these two lines for usual formulation
         grkerni = grkerni*gradhi
         grkernj = grkernj*gradh(j)
         !--use these two lines for number density formulation
         !grkerni = grkerni*(1. + gradhni*gradhi/pmassj)
         !grkernj = grkernj*(1. + gradhn(j)*gradh(j)/pmassi)
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
    prj = pr(j)        
    Prho2j = pr(j)*rho21j
    spsoundj = spsound(j)
    
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
       
!       if (geom(1:6).ne.'cartes') then
          call metric_diag(x(:,j),gdiagj(:),sqrtgj,ndim,ndimV,geom)
          B2j = dot_product_gr(Bj,Bj,gdiagj)
!       else
!          B2j = dot_product(Bj,Bj)
!       endif
       Brho2j = B2j*rho21j
       valfven2j = B2j/dens(j)
       projBconst = dot_product(Bconst,dr)
    endif
        
    !--maximum velocity for timestep control
    !            vmag = SQRT(DOT_PRODUCT(dvel,dvel))
    !            vsigdtc = vmag + spsoundi + spsoundj        &
    !                         + sqrt(valfven2i) + sqrt(valfven2j)
    vsig = 0.
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
    vsig2i = spsoundi**2 + valfven2i
    vsig2j = spsoundj**2 + valfven2j                                    
    vsigproji = vsig2i**2 - 4.*(spsoundi*(projBi))**2*rho1i
    vsigprojj = vsig2j**2 - 4.*(spsoundj*(projBj))**2*rho1j
    if (vsigproji.lt.0.) then
       write(iprint,*) ' rates: i=',i,' vsig det < 0 ',vsigproji,vsig2i**2,4*(spsoundi*projBi)**2*rho1i
       call quit 
    elseif (vsigprojj.lt.0.) then
       write(iprint,*) ' rates: j=',j,' vsig det < 0 ',vsigprojj,vsig2j**2,4*(spsoundj*projBj)**2*rho1j
       call quit  
    endif
    vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
    vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))

    vsig = vsigi + vsigj + beta*abs(dvdotr) ! also used where dvdotr>0 in MHD
    !!vsignonlin = vsigi + vsigj !!!beta*abs(dvdotr)
    
    ! vsigdtc is the signal velocity used in the timestep control
    vsigdtc = 0.5*(vsigi + vsigj + beta*abs(dvdotr))
    dvsigdtc = 1./vsigdtc
    vsigmax = max(vsigmax,vsigdtc)
    !
    !--time step control (courant and viscous)
    !
    if (vsigdtc.gt.zero) dtcourant = min(dtcourant,min(hi*dvsigdtc,hj*dvsigdtc))
    
!----------------------------------------------------------------------------
!  artificial dissipation terms
!----------------------------------------------------------------------------

    if (iav.gt.0) call artificial_dissipation

!----------------------------------------------------------------------------
!  pressure term (generalised form)
!----------------------------------------------------------------------------
    if (iprterm.ge.0) then
       if (iav.lt.0) then
          !!if (abs(pri-prj).gt.1.e-2) then
             call riemannsolver(gamma,pri,prj,-projvi,-projvj, &
               rhoi,rhoj,prstar,vstar)
!             if (prstar.gt.1.00001*max(pri,prj) .or. prstar.lt.0.99999*min(pri,prj)) then
!                print*,'error, pstar = ',prstar,'i,j = ',pri,prj,projvi,projvj
!             endif
             !prstar = 0.5*(pri + prj)
             prterm = prstar*(phii_on_phij*rho21i*sqrtgi*grkerni &
                         + phij_on_phii*rho21j*sqrtgj*grkernj)
          !!else
          !!  prterm = phii_on_phij*Prho2i*sqrtgi*grkerni &
          !!       + phij_on_phii*Prho2j*sqrtgj*grkernj
          !!endif
       else
          prterm = phii_on_phij*Prho2i*sqrtgi*grkerni &
                 + phij_on_phii*Prho2j*sqrtgj*grkernj
          !!prterm = (pri - prj)/rhoij*grkern
       endif
       !
       !--add pressure terms to force
       !
       force(:,i) = force(:,i) - pmassj*prterm*dr(:)
       force(:,j) = force(:,j) + pmassi*prterm*dr(:)
    endif

!    force(:,i) = force(:,i) + pmassj*(pr(i)/rho(i)*dri(:)*grkerni &
!                            + pr(j)/rho(j)*drj(:)*grkernj)
!    force(:,j) = force(:,j) - pmassi*(pr(i)/rho(i)*dri(:)*grkerni &
!                            + pr(j)/rho(j)*drj(:)*grkernj)

!------------------------------------------------------------------------
!  Self-gravity term
!------------------------------------------------------------------------
     
!------------------------------------------------------------------------
!  Lorentz force and time derivative of B terms
!------------------------------------------------------------------------

    if (imhd.ne.0) call mhd_terms

!------------------------------------------------------------------------
!  total energy equation (thermal energy equation terms calculated
!                         outside loop and in artificial_dissipation)
!------------------------------------------------------------------------
   
    if (iener.eq.3) then
       call energy_equation
    elseif (iener.gt.0 .and. iav.lt.0) then
       dudt(i) = dudt(i) + pmassj*prstar*(rho21i)*dvdotr*grkerni
       dudt(j) = dudt(j) + pmassi*prstar*(rho21j)*dvdotr*grkernj
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
       if (iavlim(1).eq.2) then
          !!graddivvterm = 
          graddivv(:,i) = graddivv(:,i) + pmassj*rho1j/rij*dvdotr*grkerni*dr(:)
          graddivv(:,j) = graddivv(:,j) - pmassi*rho1i/rij*dvdotr*grkernj*dr(:)
       endif
    endif

!----------------------------------------------------------------------------
!  XSPH term for moving the particles
!----------------------------------------------------------------------------
    if (ixsph.eq.1) then
       xsphterm(:,i) = xsphterm(:,i) - pmassj*dvel(:)*rhoav1*wabalt
       xsphterm(:,j) = xsphterm(:,j) + pmassi*dvel(:)*rhoav1*wabalt
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
    real :: termnonlin
    !
    !--definitions
    !      
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    alphaB = 0.5*(alphaBi + alpha(3,j))
    !!rhoav1 = 2./(rhoi + rhoj)
    if (geom(1:4).ne.'cart') then
       dpmomdotr = abs(dot_product(pmom(:,i)-pmom(:,j),dr(:)))
    else
       dpmomdotr = -dvdotr
    endif
    term = 0.5*vsig*rhoav1*grkern
    termnonlin = term
    !!termnonlin = 0.5*vsignonlin*rhoav1*grkern

    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------
    
    if (dvdotr.lt.0 .and. iav.le.2) then            
       visc = alphaav*term*dpmomdotr     ! viss=abs(dvdotr) defined in rates
       force(:,i) = force(:,i) - pmassj*visc*dr(:)
       force(:,j) = force(:,j) + pmassi*visc*dr(:)
    elseif (iav.ge.3) then ! using total energy, for approaching and receding
       visc = alphaav*term
       force(:,i) = force(:,i) + pmassj*visc*dvel(:)
       force(:,j) = force(:,j) - pmassi*visc*dvel(:)
    endif
    
    !--------------------------------------------
    !  resistivity in induction equation
    !  (applied everywhere)
    !--------------------------------------------
    if (imhd.ne.0) then
       if (iav.ge.2) then
          Bvisc(:) = dB(:)*rhoav1
       else
          Bvisc(:) = (dB(:) - dr(:)*projdB)*rhoav1 
       endif
       if (imhd.gt.0) then
          dBdtvisc(:) = alphaB*termnonlin*Bvisc(:)
       else !--vector potential resistivity
          dBdtvisc(:) = alphaB*termnonlin*(Bevol(:,i) - Bevol(:,j))
       endif
       !
       !--add to d(B/rho)/dt (converted to dB/dt later if required)
       !
       dBevoldt(:,i) = dBevoldt(:,i) + rhoi*pmassj*dBdtvisc(:)               
       dBevoldt(:,j) = dBevoldt(:,j) - rhoj*pmassi*dBdtvisc(:)
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
       if (imhd.gt.0) then
          if (iav.ge.2) then
             B2i = dot_product(Bi,Bi) ! total magnetic energy 
             B2j = dot_product(Bj,Bj) 
          else
             B2i = (dot_product(Bi,Bi) - dot_product(Bi,dr)**2) ! magnetic energy 
             B2j = (dot_product(Bj,Bj) - dot_product(Bj,dr)**2) ! along line of sight
          endif
          qdiff = qdiff + alphaB*0.5*(B2i-B2j)*rhoav1
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
          vissv = -alphaav*0.5*(dot_product(veli,dr) - dot_product(velj,dr))**2
       elseif (iav.ge.3) then
          vissv = -alphaav*0.5*dot_product(dvel,dvel)       
       else
          vissv = 0.
       endif
       !
       !  thermal energy terms
       !
       vissu = alphau*(uu(i) - uu(j))        
       !
       !  add magnetic energy term - applied everywhere
       !
       if (imhd.gt.0) then
          if (iav.ge.2) then
             vissB = -alphaB*0.5*(dot_product(dB,dB))*rhoav1 
          else
             vissB = -alphaB*0.5*(dot_product(dB,dB)-projdB**2)*rhoav1
          endif
       else
          vissB = 0.
       endif
       !
       !  add to thermal energy equation
       ! 
       dudt(i) = dudt(i) + pmassj*(term*(vissv) + termnonlin*(vissB + vissu))
       dudt(j) = dudt(j) + pmassi*(term*(vissv) + termnonlin*(vissB - vissu))    
    endif
    
    return
  end subroutine artificial_dissipation

!----------------------------------------------------------------
! Magnetic field terms - ie force and magnetic field evolution
!
!----------------------------------------------------------------
  subroutine mhd_terms
    implicit none
    real :: dalphaterm
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
       call cross_product3D(dB,dr,curlBi)
!          curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
!          curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
!          curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
!       elseif (ndimV.eq.2) then  ! just Jz in 2D
!          curlBi(1) = dB(1)*dr(2) - dB(2)*dr(1)
!          curlBi(2) = 0.
!       endif
       curlB(:,i) = curlB(:,i) + pmassj*curlBi(:)*grkern
       curlB(:,j) = curlB(:,j) + pmassi*curlBi(:)*grkern
       !
       !--add Lorentz force to total force
       !              
       select case(imagforce)
       case(1)
          fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)/rho2i
          fmag(:,j) = fmag(:,j) + pmassi*fmagj(:)/rho2j
          force(:,i) = force(:,i) + pmassj*fmagi(:)/rho2i
          force(:,j) = force(:,j) + pmassi*fmagj(:)/rho2j
       case(5)    ! Morris' Hybrid force
          fmag(:,i) = fmag(:,i) + pmassj*(faniso(:)-fiso*dr(:))
          fmag(:,j) = fmag(:,j) + pmassi*(faniso(:)+fiso*dr(:))
          force(:,i) = force(:,i) + pmassj*(faniso(:)-fiso*dr(:))
          force(:,j) = force(:,j) + pmassi*(faniso(:)+fiso*dr(:))                           
       case default      ! symmetric forces fmagxi = -fmagxj
          fmag(:,i) = fmag(:,i) + pmassj*(fmagi(:))
          fmag(:,j) = fmag(:,j) - pmassi*(fmagi(:))
          force(:,i) = force(:,i) + pmassj*(fmagi(:))
          force(:,j) = force(:,j) - pmassi*(fmagi(:))
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
       dBevoldt(:,i) = dBevoldt(:,i)            &
          + pmassj*(veli(:)*projBj + velj(:)*projBi)*rho1j*grkerni 
       dBevoldt(:,j) = dBevoldt(:,j)             &
          - pmassi*(veli(:)*projBj + velj(:)*projBi)*rho1i*grkernj
    case(3,13)   ! goes with imagforce = 3
       dBevoldt(:,i) = dBevoldt(:,i)         &
          - pmassj*projBrhoj*dvel(:)*grkerni
       dBevoldt(:,j) = dBevoldt(:,j)         &
          - pmassi*projBrhoi*dvel(:)*grkernj
    case(2,12)   ! conservative form (no change for 1D)
       dBevoldt(:,i) = dBevoldt(:,i)            &
          - phii_on_phij*pmassj*(dvel(:)*projBrhoi - rho1i*veli(:)*projdB)*grkerni 
       dBevoldt(:,j) = dBevoldt(:,j)             &
          - phij_on_phii*pmassi*(dvel(:)*projBrhoj - rho1j*velj(:)*projdB)*grkernj
    case(1,11) ! surface flux-conservative (usual) form
       dBevoldt(:,i) = dBevoldt(:,i)            &
          - phii_on_phij*pmassj*(dvel(:)*projBrhoi)*grkerni 
       dBevoldt(:,j) = dBevoldt(:,j)             &
          - phij_on_phii*pmassi*(dvel(:)*projBrhoj)*grkernj
    end select
       
    if (idivBzero.ge.2) then ! add hyperbolic correction term
       gradpsiterm = (psi(i)-psi(j))*grkern ! (-ve grad psi)
       gradpsi(:,i) = gradpsi(:,i) + pmassj*gradpsiterm*dr(:)
       gradpsi(:,j) = gradpsi(:,j) + pmassi*gradpsiterm*dr(:)
    endif

    return
  end subroutine mhd_terms
 
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
