!!--------------------------------------------------------------------
!! Computes the rates of change of the conserved variables
!! (forces, energy etc)
!! This is the core of the SPH algorithm
!!--------------------------------------------------------------------

subroutine get_rates
! USE dimen_mhd
 use debug
 use loguns
 use artvi
 use eos
 use gravity
 use hterms
 use kernel
 use linklist
 use options
 use part
 use rates
 use timestep
 use xsph
 use anticlumping
 use setup_params      ! for hfact
  
 use fmagarray
 use derivB
!
!--define local variables
!
 implicit none
 integer :: i,j,n
 integer :: icell,iprev,nneigh
 integer, allocatable, dimension(:) :: listneigh
 integer :: idone
 integer, dimension(3**ndim) :: neighcell
!
!  (particle properties - local copies and composites)
!
 real :: rij,rij2
 real :: rhoi,rho1i,rho2i,rho21i,rhoj,rho1j,rho2j,rho21j,rhoav,rhoav1,rhoij
 real :: pmassi,pmassj
 real :: Prho2i,Prho2j,prterm
 real :: hi,hi1,hj,hj1,hi21,hj21    
 real :: hav,hav1,h21
 real :: hfacwab,hfacwabi,hfacwabj,hfacgrkern,hfacgrkerni,hfacgrkernj
 real, dimension(ndim) :: dx
!
!--gr terms
!
 real :: sqrtgi,sqrtgj
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
 real :: fiso
 real :: valfven2i,valfven2j
 real :: BidotdB,BjdotdB,Brho2i,Brho2j
 real :: projBrhoi,projBrhoj,projBi,projBj,projdB
!
!  (artificial viscosity quantities)
!      
 real :: vsig,vsigi,vsigj,viss
 real :: spsoundi,spsoundj,alphai,alphaui,alphaBi
 real :: rhoi5,rhoj5
 real :: vsig2i,vsig2j,vsigproji,vsigprojj
 real :: vsigii,vsigjj
!
!  (av switch)
!
 real :: source,tdecay1,sourcedivB,sourceJ,sourceB,sourceu
 real :: graduterm
 real, dimension(:,:), allocatable :: gradu
!
!  (alternative forms)
!
 real, dimension(:), allocatable :: phi
 real :: phii,phii1,phii_on_phij,phij_on_phii
!
!  (kernel related quantities)
!
 real :: q2,q2i,q2j
 real :: wab,wabi,wabj
 real :: grkern,grkerni,grkernj
!
!  (joe's mhd fix)
!
 real :: wabaniso,wabanisoi,wabanisoj,grkernaniso,grkernanisoi,grkernanisoj
!
!  (time step criteria)
!      
 real :: vsigdtc,zero,fhmax, fonh, forcemag
!
!  (variable smoothing length terms)
!
 real :: gradhi,gradhj,gradhanisoi,gradhanisoj
 integer :: ierr
!
!  (gravity)
! 
 real :: poten
!
!--div B correction
! 
 real :: gradpsiterm,vsig2,dtcourant2

!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine get_rates'
!
!--allocate memory for local arrays
!
 nlistdim = ntotal
 allocate ( listneigh(nlistdim),STAT=ierr )
 if (ierr.ne.0) write(iprint,*) ' Error allocating neighbour list, ierr = ',ierr
 allocate ( phi(ntotal), gradu(ndim,ntotal), STAT=ierr )
 if (ierr.ne.0) write(iprint,*) ' Error allocating phi, ierr = ',ierr  
 listneigh = 0
!
!--initialise quantities
!      
 dtcourant = 1.e6  
 zero = 1.e-10
 do i=1,ntotal      ! using ntotal just makes sure they are zero for ghosts
  force(:,i) = 0.0
  drhodt(i) = 0.0
  dudt(i) = 0.0
  dendt(i) = 0.0
  dBconsdt(:,i) = 0.0
  daldt(:,i) = 0.0
  dpsidt(i) = 0.0
  fmag(:,i) = 0.0
  divB(i) = 0.0
  curlB(:,i) = 0.0
  xsphterm(:,i) = 0.0
  gradu(:,i) = 0.0
 enddo
!
!--set MHD quantities to zero if mhd not set
!
 if (imhd.le.0) then  ! these quantities are still used if mhd off
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
 endif

!
!--skip the whole neighbour thing if it is doing nothing
!
 if (iprterm.lt.0 .and. iav.eq.0 .and. imhd.eq.0 .and. iener.eq.0 &
      .and. ihvar.lt.2 .and. icty.eq.0) then
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
!
!--Loop over all the link-list cells
!
 loop_over_cells: do icell=1,ncellsloop            ! step through all cells

    !print*,'> doing cell ',icell
    !
    !--get the list of neighbours for this cell 
    !  (common to all particles in the cell)
    !
    call get_neighbour_list(icell,neighcell,listneigh,nneigh)
       
    i = ifirstincell(icell)      ! start with first particle in cell
    idone = -1
    if (i.ne.-1) iprev = i

    loop_over_cell_particles: do while (i.ne.-1)      ! loop over home cell particles

       !print*,'Doing particle ',i,x(:,i),valfven2i,rho(i),hh(i)
       idone = idone + 1
       rhoi = rho(i)
       rho2i = rhoi*rhoi
       rhoi5 = sqrt(rhoi)
       rho1i = 1./rhoi
       rho21i = rho1i*rho1i       
       Prho2i = pr(i)*rho21i
       spsoundi = spsound(i)
       veli(:) = vel(:,i)
       pmassi = pmass(i)
       alphai = alpha(1,i)
       alphaui = alpha(2,i)
       phii = phi(i)
       phii1 = 1./phii
       sqrtgi = sqrtg(i)
       if (imhd.ge.11) then      ! if mag field variable is B
          Bi(:) = Bcons(:,i)
          Brhoi(:) = Bi(:)*rho1i
       elseif (imhd.gt.0) then      ! if mag field variable is B/rho
          Brhoi(:) = Bcons(:,i)
          Bi(:) = Bfield(:,i)
       endif
       ! mhd definitions
       if (imhd.ne.0) then
          Brho2i = dot_product(Brhoi,Brhoi)
          valfven2i = Brho2i*rhoi
          alphaBi = alpha(3,i)
       endif
       gradhi = 1./(1. - gradh(i))
       gradhanisoi = 1./(1.-gradhaniso(i))
       !if (gradhi.le.0.5) then
       !   write(iprint,*) 'Error in grad h terms, part ',i,gradhi
       !endif
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
           hfacwabj = hj1**ndim
           hfacgrkernj = hfacwabj*hj1
           hav = 0.5*(hi + hj)
           hav1 = 1./hav
           h21 = hav1*hav1
           hfacwab = hav1**ndim 
           hfacgrkern = hfacwab*hav1
     
           rij2 = dot_product(dx,dx)
           rij = sqrt(rij2)
           if (rij.le.1.e-10) then
              write(iprint,*) 'rates: dx = 0 i,j,dx,hi,hj=',i,j,dx,hi,hj
                call quit
           endif      
           q2 = rij2*h21      
           q2i = rij2*hi21
           q2j = rij2*hj21
           dr(1:ndim) = dx(1:ndim)/rij      ! unit vector
           if (ndimV.gt.ndim) dr(ndim+1:ndimV) = 0.

             !----------------------------------------------------------------------------
             !  do pairwise interaction if either particle is within range of the other
             !----------------------------------------------------------------------------
           if ((q2i.lt.radkern2).or.(q2j.lt.radkern2)) then  ! if < 2h
                call rates_core
           else      ! if outside 2h
!              PRINT*,'outside 2h, not calculated, r/h=',sqrt(q2)
           endif      ! r < 2h or not

         endif            ! j.ne.i
        
        enddo loop_over_neighbours

       iprev = i
       if (iprev.ne.-1) i = ll(i)            ! possibly should be only IF (iprev.NE.-1)
    
    enddo loop_over_cell_particles
            
 enddo loop_over_cells
 
666 continue

!----------------------------------------------------------------------------
!  calculate gravitational force on all the particles
!----------------------------------------------------------------------------
 if (igravity.ne.0) call direct_sum_poisson( &
                     x(:,1:npart),pmass(1:npart),poten,fgrav(:,1:npart),npart)

 if (trace) write(iprint,*) 'Finished main rates loop'
 fhmax = 0.0
 dtforce = 1.e6
 dtcourant2 = 1.e6

!----------------------------------------------------------------------------
!  loop over the particles again, subtracting external forces and source terms
!----------------------------------------------------------------------------
 
 do i=1,npart
 
    rho1i = 1./rho(i)
!
!--subtract external forces
!
    if (iexternal_force.ne.0) call external_forces(i,iexternal_force)
!
!--add self-gravity force
!
    if (igravity.ne.0) force(1:ndim,i) = force(1:ndim,i) - fgrav(1:ndim,i)*rho1i
!
!--damp force if appropriate
!
    if (damp.gt.0.) force(:,i) = force(:,i) - damp*vel(:,i)
!
!--add source terms (derivatives of metric) to momentum equation
!
    if (igeom.ne.0 .and. allocated(sourceterms)) force(:,i) = force(:,i) + sourceterms(:,i)
!
!--do the divisions by rho etc (this is for speed - so calculations are not
!  done multiple times within the loop)
!
    if (imhd.ne.0) then
       divB(i) = divB(i)*rho1i            !*rhoi
       curlB(:,i) = curlB(:,i)*rho1i
       dBconsdt(:,i) = dBconsdt(:,i)*rho1i
    endif    
    gradu(:,i) = gradu(:,i)*rho1i
!
!--calculate time derivative of the smoothing length from the density derivative
!
    if (ihvar.eq.2 .or. ihvar.eq.3) then
       dhdt(i) = -hh(i)/(ndim*rho(i))*drhodt(i)
    else
       dhdt(i) = 0.    
    endif
!
!--if using the thermal energy equation, set the energy derivative
!  (note that dissipative terms are calculated in rates, but otherwise comes straight from cty)
!
    if (iener.ne.3 .and. iener.ne.0) then
       dudt(i) = dudt(i) + pr(i)*rho1i**2*drhodt(i)
       dendt(i) = dudt(i)
    endif
!
!--calculate maximum force/h for the timestep condition
!  also check for errors in the force
!
    if ( any(force(:,i).gt.1.e8)) then
       write(iprint,*) 'rates: force ridiculous ',force(:,i)
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
    if (imhd.ne.0) valfven2i = dot_product(Bfield(:,i),Bfield(:,i))*rho1i
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
       if (iavlim(1).ne.0) then
          source = max(drhodt(i)*rho1i,0.0)
          if (iavlim(1).eq.2) source = source*(2.0-alpha(1,i))      
          daldt(1,i) = (alphamin - alpha(1,i))*tdecay1 + avfact*source
       endif
       !
       !--artificial thermal conductivity parameter
       !
       if (iavlim(2).ne.0 .and. iener.gt.0) then
          if (uu(i).gt.0.) then
	     sourceu = 0.5*sqrt(dot_product(gradu(:,i),gradu(:,i))/uu(i))
          else
	     sourceu = 0.
	  endif
	  if (iavlim(2).eq.2) sourceu = sourceu*(2.0-alpha(2,i)) 
          daldt(2,i) = (alphaumin - alpha(2,i))*tdecay1 + sourceu       
       endif
       !
       !--artificial resistivity parameter if iavlim > 10
       !
       if (iavlim(3).ne.0 .and. imhd.gt.0) then
          !--calculate source term for the resistivity parameter
          sourceJ = SQRT(DOT_PRODUCT(curlB(:,i),curlB(:,i))*rho1i)
          sourcedivB = abs(divB(i))*SQRT(rho1i)
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
          dpsidt(i) = -vsig2*divB(i) - psidecayfact*psi(i)*vsig/hh(i)          
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
       dBconsdt(:,i) = 0.0
       daldt(:,i) = 0.
       dpsidt(i) = 0.0
       fmag(:,i) = 0.0
       divB(i) = 0.0
       curlB(:,i) = 0.0
       xsphterm(:,i) = 0.0
    endif
 enddo

 if (allocated(listneigh)) deallocate(listneigh)
 if (allocated(phi)) deallocate(phi,gradu)
 if (trace) write(iprint,*) ' Exiting subroutine get_rates'
      
 return

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
!!    use unityfunc
    implicit none
!
!--calculate both kernels if using anticlumping kernel
!
   if (imhd.ne.0 .and. imagforce.eq.2 .and. ianticlump.eq.1) then
      if (ikernav.eq.1) then
         call interpolate_kernels(q2,wab,grkern,wabaniso,grkernaniso)
         wab = wab*hfacwab
         grkern = grkern*hfacgrkern
         grkernaniso = grkernaniso*hfacgrkern
         grkerni = grkern
         grkernj = grkern
         grkernanisoi = grkernaniso
         grkernanisoj = grkernaniso
      else
         !  (using hi)
         call interpolate_kernels(q2i,wabi,grkerni,wabanisoi,grkernanisoi)
         wabi = wabi*hfacwabi
         grkerni = grkerni*hfacgrkerni
         grkernanisoi = grkernanisoi*hfacgrkerni
         !  (using hj)
         call interpolate_kernels(q2j,wabj,grkernj,wabanisoj,grkernanisoj)
         wabj = wabj*hfacwabj
         grkernj = grkernj*hfacgrkernj
         grkernanisoj = grkernanisoj*hfacgrkernj
         !  (calculate average)              
         wab = 0.5*(wabi + wabj)  ! wab only used in XSPH term
         !  (grad h terms)  
         if (ikernav.eq.3) then  ! if using grad h correction
            gradhj = 1./(1. - gradh(j))
            grkerni = grkerni*gradhi
            grkernj = grkernj*gradhj
            grkern = 0.5*(grkerni + grkernj)
            gradhanisoj = 1./(1. - gradhaniso(j))
            grkernanisoi = grkernanisoi*gradhanisoi
            grkernanisoj = grkernanisoj*gradhanisoj
            !!grkernaniso = 0.5*(grkernanisoi + grkernanisoj)	  
         else  ! if not using grad h correction               
            grkern = 0.5*(grkerni + grkernj)
            grkerni = grkern
            grkernj = grkern
            grkernaniso = 0.5*(grkernanisoi + grkernanisoj)
            grkernanisoi = grkernaniso
            grkernanisoj = grkernaniso
         endif
      endif
   else
!
!--or just calculate the usual kernel
!
      if (ikernav.eq.1) then
         call interpolate_kernel(q2,wab,grkern)
         wab = wab*hfacwab
         grkern = grkern*hfacgrkern
         grkerni = grkern
         grkernj = grkern
      else
         !  (using hi)
         call interpolate_kernel(q2i,wabi,grkerni)
         wabi = wabi*hfacwabi
         grkerni = grkerni*hfacgrkerni
         !  (using hj)
         call interpolate_kernel(q2j,wabj,grkernj)
         wabj = wabj*hfacwabj
         grkernj = grkernj*hfacgrkernj
         !  (calculate average)              
         wab = 0.5*(wabi + wabj)
         !  (grad h terms)  
         if (ikernav.eq.3) then  ! if using grad h correction
            gradhj = 1./(1. - gradh(j))
            grkerni = grkerni*gradhi
            grkernj = grkernj*gradhj
            grkern = 0.5*(grkerni + grkernj)
         else  ! if not using grad h correction               
            grkern = 0.5*(grkerni + grkernj)
            grkerni = grkern
            grkernj = grkern
         endif
      endif
   endif
!
!--define local copies of quantities
!
    velj(:) = vel(:,j)            
    dvel(:) = veli(:) - velj(:)
    dvdotr = dot_product(dvel,dr)
    rhoj = rho(j)
    rho1j = 1./rhoj
    rho2j = rhoj*rhoj
    rho21j = rho1j*rho1j
    rhoj5 = sqrt(rhoj)
    rhoij = rhoi*rhoj            
    Prho2j = pr(j)*rho21j
    spsoundj = spsound(j)
    pmassj = pmass(j)
    
    phii_on_phij = phii/phi(j)
    phij_on_phii = phi(j)*phii1       
    sqrtgj = sqrtg(j)      
    !-- mhd definitions --!
    if (imhd.ge.11) then      ! if B is mag field variable
       Bj(:) = Bcons(:,j)
       Brhoj(:) = Bj(:)*rho1j
    elseif (imhd.ne.0) then      ! if B/rho is mag field variable
       Brhoj(:) = Bcons(:,j)
       Bj(:) = Bfield(:,j)                                
    endif
    if (imhd.ne.0) then
       dB(:) = Bi(:) - Bj(:)
       projBi = dot_product(Bi,dr)
       projBj = dot_product(Bj,dr)
       projdB = dot_product(dB,dr)
       projBrhoi = dot_product(Brhoi,dr)
       projBrhoj = dot_product(Brhoj,dr)
       Brho2j = dot_product(Brhoj,Brhoj)
       valfven2j = Brho2j*rhoj
    endif
        
    !--maximum velocity for timestep control
    !            vmag = SQRT(DOT_PRODUCT(dvel,dvel))
    !            vsigdtc = vmag + spsoundi + spsoundj        &
    !                         + sqrt(valfven2i) + sqrt(valfven2j)
    vsig = 0.
    viss = 0.
!----------------------------------------------------------------------------
!  calculate signal velocity (this is used for timestep control
!  and also in the artificial viscosity)
!----------------------------------------------------------------------------

    if (dvdotr.lt.0) viss = abs(dvdotr)
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
    vsigproji = vsig2i**2 - 4.*(spsoundi*projBi)**2*rho1i
    vsigprojj = vsig2j**2 - 4.*(spsoundj*projBj)**2*rho1j
    if (vsigproji.lt.0. .or. vsigprojj.lt.0.) then
       write(iprint,*) ' rates: vsig det < 0 ', &
       'i: ',vsigproji,vsig2i**2,4*(spsoundi*projBi)**2*rho1i,      &
       'j: ',vsigprojj,vsig2j**2,4*(spsoundj*projBj)**2*rho1j
       call quit  
    endif
    vsigi = SQRT(0.5*(vsig2i + SQRT(vsigproji)))
    vsigj = SQRT(0.5*(vsig2j + SQRT(vsigprojj)))

    vsig = vsigi + vsigj + beta*abs(dvdotr) ! also used where dvdotr>0 in MHD
    
    ! vsigdtc is the signal velocity used in the timestep control
    vsigdtc = vsigi + vsigj + beta*abs(dvdotr)
    !
    !--time step control (courant and viscous)
    !
    if (vsigdtc.gt.zero) dtcourant = min(dtcourant,min(hi,hj)/vsigdtc)
    
!----------------------------------------------------------------------------
!  artificial dissipation terms
!----------------------------------------------------------------------------

    if (iav.ne.0) call artificial_dissipation

!----------------------------------------------------------------------------
!  pressure term (generalised form)
!----------------------------------------------------------------------------
    if (iprterm.ge.0) then
       prterm = phii_on_phij*Prho2i*sqrtgi*grkerni &
              + phij_on_phii*Prho2j*sqrtgj*grkernj
    else
       prterm = 0.
    endif
    !
    !--add pressure terms to force
    !
    force(:,i) = force(:,i) - pmassj*prterm*dr(:) !!!- grkerncorri*gradunity(:,i))
    force(:,j) = force(:,j) + pmassi*prterm*dr(:) !!!- grkerncorrj*gradunity(:,j))
!----------------------------------------------------------------------------
!  time derivative of density (continuity equation) in generalised form
!  compute this even if direct sum - used in dudt, dhdt and av switch
!----------------------------------------------------------------------------
!
    !drhodt(i) = drhodt(i) + phii_on_phij*pmassj*dvdotr*grkerni
    !drhodt(j) = drhodt(j) + phij_on_phii*pmassi*dvdotr*grkernj    
    drhodt(i) = drhodt(i) + pmassj*dvdotr*grkerni
    drhodt(j) = drhodt(j) + pmassi*dvdotr*grkernj

!------------------------------------------------------------------------
!  Lorentz force and time derivative of B terms
!------------------------------------------------------------------------

    if (imhd.ne.0) call mhd_terms

!------------------------------------------------------------------------
!  total energy equation (thermal energy equation terms calculated
!                         outside loop and in artificial_dissipation)
!------------------------------------------------------------------------
   
    if (iener.eq.3) call energy_equation

!------------------------------------------------------------------------
!  grad u term for dissipation switch
!------------------------------------------------------------------------

    if (iav.ne.0 .and. iavlim(2).ne.0) then
       graduterm = uu(j)-uu(i)
       gradu(:,i) = gradu(:,i) + pmassj*graduterm*grkerni*dr(:)
       gradu(:,j) = gradu(:,j) + pmassi*graduterm*grkernj*dr(:)
    endif

!----------------------------------------------------------------------------
!  XSPH term for moving the particles
!----------------------------------------------------------------------------
    if (ixsph.eq.1) then
       xsphterm(:,i) = xsphterm(:,i) - pmassj*dvel(:)*rhoav1*wab
       xsphterm(:,j) = xsphterm(:,j) + pmassi*dvel(:)*rhoav1*wab
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
    real :: term,rhoav1
    !
    !--definitions
    !      
    alphaav = 0.5*(alphai + alpha(1,j))
    alphau = 0.5*(alphaui + alpha(2,j))
    alphaB = 0.5*(alphaBi + alpha(3,j))
    rhoav1 = 2./(rhoi + rhoj)
    term = 0.5*vsig*rhoav1*grkern

    !----------------------------------------------------------------
    !  artificial viscosity in force equation
    ! (applied only for approaching particles)
    !----------------------------------------------------------------
    
    if (dvdotr.lt.0 .and. iav.le.2) then            
       visc = alphaav*term*viss     ! viss=abs(dvdotr) defined in rates
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
    if (iav.ge.2) then
       Bvisc(:) = dB(:)*rhoav1
    else
       Bvisc(:) = (dB(:) - dr(:)*projdB)*rhoav1 
    endif
    dBdtvisc(:) = alphaB*term*Bvisc(:)
    
    if (imhd.le.10) then            ! evolving B/rho
       dBconsdt(:,i) = dBconsdt(:,i) + rhoi*pmassj*dBdtvisc(:)               
       dBconsdt(:,j) = dBconsdt(:,j) - rhoj*pmassi*dBdtvisc(:)
    elseif (imhd.eq.11) then      ! evolving B
       dBconsdt(:,i) = dBconsdt(:,i) + rho2i*pmassj*dBdtvisc(:)               
       dBconsdt(:,j) = dBconsdt(:,j) - rho2j*pmassi*dBdtvisc(:)
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
       if (iav.ge.2) then
          B2i = dot_product(Bi,Bi) ! total magnetic energy 
          B2j = dot_product(Bj,Bj) 
       else
          B2i = (dot_product(Bi,Bi) - dot_product(Bi,dr)**2) ! magnetic energy 
          B2j = (dot_product(Bj,Bj) - dot_product(Bj,dr)**2) ! along line of sight
       endif
       qdiff = qdiff + alphaB*0.5*(B2i-B2j)*rhoav1
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
       if (imhd.ne.0) then
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
       dudt(i) = dudt(i) + pmassj*term*(vissv + vissB + vissu)
       dudt(j) = dudt(j) + pmassi*term*(vissv + vissB - vissu)    
    endif
    
    return
  end subroutine artificial_dissipation

!----------------------------------------------------------------
! Magnetic field terms - ie force and magnetic field evolution
!
!----------------------------------------------------------------
  subroutine mhd_terms
    implicit none
    real :: divBonrho
    !----------------------------------------------------------------------------            
    !  Lorentz force
    !----------------------------------------------------------------------------
    !
    !--calculate curl B for current (only if field is 3D)
    !  this is used in the switch for the artificial resistivity term
    !
    if (ndimB.eq.3) then
       curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
       curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
       curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
    elseif (ndimB.eq.2) then  ! just Jz in 2D
       curlBi(1) = dB(1)*dr(2) - dB(2)*dr(1)
       curlBi(2) = 0.
    endif
    
    if (imagforce.eq.1) then      ! vector form (dot products)
       
       BidotdB = dot_product(Brhoi,dB)
       BjdotdB = dot_product(Brhoj,dB)
       fmagi(:) = grkern*(dr(:)*BidotdB - projBrhoi*dB(:))              
       fmagj(:) = grkern*(dr(:)*BjdotdB - projBrhoj*dB(:))
       
    elseif (imagforce.eq.2) then      ! tensor formalism in generalised form
       !
       !--isotropic mag force (pressure gradient)
       !              
       fiso = 0.5*(Brho2i*phij_on_phii*grkerni + Brho2j*phii_on_phij*grkernj)
       !
       !--anisotropic force (including variable smoothing length terms)
       !  (also including Joe's correction term which preserves momentum conservation)
       !   in 1D we only do this in the x-direction
       !
       if (ianticlump.eq.1) then
          faniso(1:ndim) = Brhoi(1:ndim)*projBrhoi*phij_on_phii*grkernanisoi  &
                         + Brhoj(1:ndim)*projBrhoj*phii_on_phij*grkernanisoj               
          
!          faniso(1:ndim) = Brhoi(1:ndim)*projBrhoi*phij_on_phii*grkerni  &
!                         + Brhoj(1:ndim)*projBrhoj*phii_on_phij*grkernj
!	  divBonrho =  projBrhoi*rho1i*grkerni + &
!	               projBrhoj*rho1j*grkernj
	  
          if (ndimV.gt.ndim) then
             faniso(ndim+1:ndimV) = Brhoi(ndim+1:ndimV)*projBrhoi*phij_on_phii*grkerni  &
                                  + Brhoj(ndim+1:ndimV)*projBrhoj*phii_on_phij*grkernj               
          endif

       else
          faniso(:) = Brhoi(:)*projBrhoi*phij_on_phii*grkerni  &
                    + Brhoj(:)*projBrhoj*phii_on_phij*grkernj               
       endif

       !
       !--add contributions to magnetic force
       !        
       fmagi(:) = faniso(:) - fiso*dr(:)
    
    elseif (imagforce.eq.3) then  ! alternative symmetric formulation on aniso

       rhoij = rhoi*rhoj
       fiso = 0.5*(Brho2i*grkerni + Brho2j*grkernj)
       ! faniso is        B*div B        +    B dot grad B
       faniso(:) = (Bi(:)*projBj*grkerni + Bj(:)*projBi*grkernj)/rhoij
       fmagi(:) = faniso(:) - fiso*dr(:)   
    
    elseif (imagforce.eq.5) then      ! Morris' Hybrid form
       
       rhoij = rhoi*rhoj
       fiso = 0.5*(Brho2i*grkerni + Brho2j*grkernj)
!--note that I have tried faniso with grkerni and grkernj but much worse on mshk2
       faniso(:) = grkern*(Bj(:)*projBj - Bi(:)*projBi)/rhoij
       fmagi(:) = faniso(:) - fiso*dr(:)              
       
    endif
    !
    !--compute rho * divergence of B
    !
    divB(i) = divB(i) - pmassj*projdB*grkern
    divB(j) = divB(j) - pmassi*projdB*grkern
    !divB(i) = divB(i) + pmassj*projBrhoj*grkerni
    !divB(j) = divB(j) - pmassi*projBrhoi*grkernj
    !
    !--compute rho * current density J
    !
    curlB(:,i) = curlB(:,i) - pmassj*curlBi(:)*grkern
    curlB(:,j) = curlB(:,j) - pmassi*curlBi(:)*grkern
    !
    !--add Lorentz force to total force
    !              
    if (imagforce.eq.1) then
       fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)/rho2i
       fmag(:,j) = fmag(:,j) - pmassi*fmagj(:)/rho2j
       force(:,i) = force(:,i) + pmassj*fmagi(:)/rho2i
       force(:,j) = force(:,j) - pmassi*fmagj(:)/rho2j
    elseif (imagforce.eq.5) then      ! Morris' Hybrid force
       fmag(:,i) = fmag(:,i) + pmassj*(faniso(:)-fiso*dr(:))
       fmag(:,j) = fmag(:,j) + pmassi*(faniso(:)+fiso*dr(:))
       force(:,i) = force(:,i) + pmassj*(faniso(:)-fiso*dr(:))
       force(:,j) = force(:,j) + pmassi*(faniso(:)+fiso*dr(:))                           
    else       ! symmetric forces fmagxi = -fmagxj
       fmag(:,i) = fmag(:,i) + pmassj*(fmagi(:))
       fmag(:,j) = fmag(:,j) - pmassi*(fmagi(:))
       force(:,i) = force(:,i) + pmassj*(fmagi(:))
       force(:,j) = force(:,j) - pmassi*(fmagi(:))
    endif
    
    !--------------------------------------------------------------------------------
    !  time derivative of magnetic field (divide by rho later) - in generalised form
    !---------------------------------------------------------------------------------
    !
    !   (evolving B/rho)
    !
    if (imhd.gt.0.and. imhd.le.10) then   ! divided by rho later
       if (imhd.eq.4) then   ! conservative form (explicitly symmetric)
          dBconsdt(:,i) = dBconsdt(:,i)            &
               + pmassj*(veli(:)*projBj + velj(:)*projBi)*rho1j*grkerni 
          dBconsdt(:,j) = dBconsdt(:,j)             &
               - pmassi*(veli(:)*projBj + velj(:)*projBi)*rho1i*grkernj
       elseif (imhd.eq.3) then   ! goes with imagforce = 3
          dBconsdt(:,i) = dBconsdt(:,i)         &
               - pmassj*projBrhoj*dvel(:)*grkerni 
          dBconsdt(:,j) = dBconsdt(:,j)         &
               - pmassi*projBrhoi*dvel(:)*grkernj            
       elseif (imhd.eq.2) then  ! conservative form (no change for 1D)
          dBconsdt(:,i) = dBconsdt(:,i)            &
               - phii_on_phij*pmassj*(dvel(:)*projBrhoi - rho1i*veli(:)*projdB)*grkerni 
          dBconsdt(:,j) = dBconsdt(:,j)             &
               - phij_on_phii*pmassi*(dvel(:)*projBrhoj - rho1j*velj(:)*projdB)*grkernj
       else                 ! non-conservative (usual) form
          dBconsdt(:,i) = dBconsdt(:,i)            &
               - phii_on_phij*pmassj*(dvel(:)*projBrhoi)*grkerni 
          dBconsdt(:,j) = dBconsdt(:,j)             &
               - phij_on_phii*pmassi*(dvel(:)*projBrhoj)*grkernj
       endif
       
       if (idivBzero.ge.2) then ! add hyperbolic correction term
          gradpsiterm = (psi(i)-psi(j))*grkern ! (-ve grad psi)
          dBconsdt(:,i) = dBconsdt(:,i) + rho1i*pmassj*gradpsiterm*dr(:)
          dBconsdt(:,j) = dBconsdt(:,j) + rho1j*pmassi*gradpsiterm*dr(:)
       endif
       !
       !   (evolving B)
       !
    elseif (imhd.ge.11) then      ! note divided by rho later              
       
       dBconsdt(:,i) = dBconsdt(:,i) + pmassj*(Bi(:)*dvdotr - dvel(:)*projBi)*grkerni   
       dBconsdt(:,j) = dBconsdt(:,j) + pmassi*(Bj(:)*dvdotr - dvel(:)*projBj)*grkernj
        
       if (imhd.eq.12) then ! add v div B term
        dBconsdt(:,i) = dBconsdt(:,i) + pmassj*(-veli(:)*projdB)*grkerni   
        dBconsdt(:,j) = dBconsdt(:,j) + pmassi*(-velj(:)*projdB)*grkernj
       endif
       
       if (idivBzero.ge.2) then  ! add hyperbolic correction term
          gradpsiterm = (psi(i)-psi(j))*grkern ! (-ve rho*grad psi)   
          dBconsdt(:,i) = dBconsdt(:,i) + pmassj*gradpsiterm*dr(:)
          dBconsdt(:,j) = dBconsdt(:,j) + pmassi*gradpsiterm*dr(:)
       endif
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
    if (imhd.ne.0 .and. ianticlump.eq.1 .and. imagforce.eq.2) then
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
