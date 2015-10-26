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

!---------------------------------------------------------------------
! these subroutines calculates the primitive variables from the
! conservative variables and vice versa.
!
! This version is for non-relativistic MHD.
! The thermal energy is calculated from the
! total energy, whilst the magnetic field B is calculated from
! the conserved variable B/rho.
!
! These subroutines would be used in the GR case, for all the 
! variables. The evaluation would be far more complicated, however.
!---------------------------------------------------------------------
module cons2prim
 implicit none
 logical, parameter :: specialrelativity = .false.
 integer, parameter, private :: Jcurltype = 2

contains

subroutine conservative2primitive
  use dimen_mhd
  use debug
  use options
  use loguns
  use bound
  use eos
  use part
  use getcurl,       only:get_curl
  use getBeulerpots, only:get_B_eulerpots
  use getlambda,     only:get_lambda
  use get2ndderivs,  only:get_2ndderivs
  use smooth, only:smooth_variable
  use rates,  only:gradpsi
  use hterms, only:gradgradh,gradh,zeta
  use derivB, only:curlB,gradB
  use timestep, only:nsteps
  use part_in, only:Bevolin,dustevolin !,velin
  !use utils,   only:curl3D_epsijk
  !use khsetup, only:densmedium,denszero,smoothl,yprofile,przero
  implicit none
  integer :: i,j,nerr,nerr1,k
  real :: B2i, v2i, pri, dhdrhoi, emag, emagold, dx
  real, dimension(ndimV) :: Binti,Bfieldi !,curlBi
  real, dimension(ndim) :: dxbound
  logical, parameter :: JincludesBext = .true.
  logical :: remap
  real, parameter :: remap_tol = 1.e-1

  if (trace) write(iprint,*) ' Entering subroutine conservative2primitive'

  nerr = 0
  nerr1 = 0
  sqrtg = 1.
  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     !--error checking on dust-to-gas ratio
     do i=1,npart
        if (use_sqrtdustfrac) then
           dustfrac(i) = dustevol(i)**2/rho(i)
        else
           dustfrac(i) = dustevol(i)
        endif
        if (dustfrac(i) < 0.) then
           nerr = nerr + 1
           dustfrac(i) = 0. !abs(dustfrac(i))
           if (.not. use_sqrtdustfrac) then
              dustevol(i) = 0.
              dustevolin(i) = 0.
           endif
           !call quit
        elseif (dustfrac(i) > 1.) then
           nerr1 = nerr1 + 1
           dustfrac(i) = 1. - epsilon(1.)
           if (.not. use_sqrtdustfrac) then
              dustevol(i) = 1. - epsilon(1.)
              dustevolin(i) = 1. - epsilon(1.)
           endif
        endif
     enddo
     dens = rho*(1. - dustfrac) ! rho = rho_gas + rho_dust, dens is rho_gas
  else
     dens = rho
  endif
  if (nerr > 0) print*,'ERROR: dust fraction < 0 on ',nerr,' particles'
  if (nerr1 > 0) print*,'ERROR: dust fraction > 1 on ',nerr1,' particles'
! 
!--set x0 correctly on ghost particles
!  (needed for remapped B or remapped euler potentials evolution)
!
  if (imhd.eq.-3 .or. imhd.eq.10 .or. imhd.eq.20 .and. any(ibound.gt.1)) then
     dxbound(:) = xmax(:) - xmin(:)
     do i=npart+1,ntotal
        j = ireal(i)
        do k=1,ndim
           dx = x(k,j) - x(k,i)
           if (abs(dx).gt.0.5*dxbound(k)) then
              x0(k,i) = x0(k,j) - dxbound(k)*SIGN(1.,dx)
           else
              x0(k,i) = x0(k,j)              
           endif
        enddo
     enddo
  endif

!
!--calculate magnetic flux density B from the conserved variable
!
  remap = .false.
  if ((imhd.eq.10 .or. imhd.eq.20 .or. imhd.eq.-3) .and. &
      nsteps_remap.gt.0 .and. mod(nsteps,nsteps_remap).eq.0) then
     remap = .true.
     !print*,' REMAPPING...'
  endif

  select case(imhd)
  case(11:19, 21:) ! if using B as conserved variable
     Bfield = Bevol
  case(20)  ! remapped B
     !--remap B_0 to current B
     call get_B_eulerpots(4,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap,remap_tol)
     if (remap) then
        !
        !--recompute B field with remapped potentials (if remap=.true. on first
        !  call then Bevol comes out changed)
        !
        emagold = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy before remapping = ',emagold

        Bfield = Bevol  ! set B field equal to remapped B_0
        emag = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy after remapping = ',emag, ' change = ',(emag-emagold)/emagold
     else
        print*,' magnetic energy = ',emag_calc(pmass,rho,Bfield,npart)
     endif
  case(10)  ! remapped B/rho
     !--remap B/rho to current B/rho
     call get_B_eulerpots(3,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap,remap_tol)
     do i=1,npart
        Bfield(:,i) = Bfield(:,i)*rho(i)
     enddo
     if (remap) then
        !
        !--recompute B field with remapped potentials (if remap=.true. on first
        !  call then Bevol comes out changed)
        !
        emagold = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy before remapping = ',emagold
        do i=1,npart
           Bfield(:,i) = Bevol(:,i)*rho(i)
        enddo
        emag = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy after remapping = ',emag, ' change = ',(emag-emagold)/emagold
     else
        print*,' magnetic energy = ',emag_calc(pmass,rho,Bfield,npart)
     endif
  case(1:9) ! if using B/rho as conserved variable
     do i=1,npart
        Bfield(:,i) = Bevol(:,i)*rho(i)
     enddo
  case(-1,-2) ! if using vector potential
     !write(iprint,*) 'getting B field from vector potential...'
     if (iprterm.ge.10) stop 'conflict with use of psi variable (iprterm=10/imhd<0)'
     gradpsi(:,:) = 0.
     psi(:) = 0.
     zeta(:) = 0.
     if (iuse_exact_derivs.gt.0) then
        !--get an exact linear curl A (iderivtype=2 input to get_B_eulerpots)
        call get_B_eulerpots(2,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap)
     else
        call get_curl(1,npart,x,pmass,rho,hh,Bevol,Bfield,gradpsi)
     endif
     do i=1,npart
        Binti(:) = Bfield(:,i)
        Bfieldi(:) = Binti(:) + Bconst(:)
        if (JincludesBext) Bfield(:,i) = Bfield(:,i) + Bconst(:)
        if (imagforce.eq.7 .and. ihvar.gt.0) then
           !--construct xi term
           dhdrhoi = -hh(i)/(ndim*(rho(i)))
           gradpsi(:,i) = dhdrhoi*gradpsi(:,i)
           zeta(i) = dot_product(Bfieldi(:),gradpsi(:,i)) + dot_product(Bfieldi(:),Binti(:))*gradgradh(i)*gradh(i)
           psi(i) = dot_product(Bfieldi(:),Binti(:))*gradh(i)/rho(i)*dhdrhoi
           !if (i.lt.10) print*,i,'xi = ',psi(i),' zeta = ',gradgradh(i),' gradpsi = ',gradpsi(:,i) !,Bfield(:,i)
        else
           zeta(i) = 0.
           psi(i) = 0.
        endif
     enddo
     !if (imagforce.eq.7) then
        !--copy Bfield onto ghosts
        if (any(ibound.eq.1)) then
           do i=1,npart
              if (itype(i).eq.itypebnd) then
                 j = ireal(i)
                 if (j > 0) Bfield(:,i) = Bfield(:,j)
              endif
           enddo
        endif
        if (any(ibound.gt.1)) then
           do i=npart+1,ntotal
              j = ireal(i)
              if (j > 0) Bfield(:,i) = Bfield(:,j)
           enddo
        endif
        !--get J using symmetric curl operator
        call get_curl(Jcurltype,npart,x,pmass,rho,hh,Bfield,curlB)
     !endif
     if (.not.JincludesBext) then
        do i=1,npart
           Bfield(:,i) = Bfield(:,i) + Bconst(:)
        enddo
     endif

     !--reset gradpsi to zero after we have finished using it
     gradpsi(:,:) = 0.
  case(:-3) ! generalised Euler potentials
     write(iprint,*) 'getting B field from Generalised Euler Potentials... '
     call get_B_eulerpots(1,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap,remap_tol)
     print*,' magnetic energy = ',emag_calc(pmass,rho,Bfield,npart)
     !--add constant field component
     do i=1,npart
        Bfield(:,i) = Bfield(:,i) + Bconst(:)
     enddo

     if (remap) then
        !
        !--recompute B field with remapped potentials
        !
        emagold = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy before remapping = ',emagold
        remap = .false.
        call get_B_eulerpots(1,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap)
        remap = .true.
        emag = emag_calc(pmass,rho,Bfield,npart)
        print*,' magnetic energy after remapping = ',emag, ' change = ',(emag-emagold)/emagold
     endif

     !--copy Bfield onto ghosts
     if (any(ibound.eq.1)) then
        do i=1,npart
           if (itype(i).eq.itypebnd) then
              j = ireal(i)
              Bfield(:,i) = Bfield(:,j)
           endif
        enddo
     endif
     if (any(ibound.gt.1)) then
        do i=npart+1,ntotal
           j = ireal(i)
           Bfield(:,i) = Bfield(:,j)
        enddo
     endif
  case default
     !--no magnetic field
  end select
!
!--if magnetic field has been remapped, reset Bevolin for the timestepping
!
  if (remap) then
     Bevolin = Bevol
  endif
!
!--get magnetic current, needed for ambipolar diffusion calculation
!
  if (iambipolar > 0 .or. iavlim(3)==2) then
     if (iavlim(3)==2) then
        !--get J and grad(B) using standard (differenced) curl operator
        call get_curl(1,npart,x,pmass,rho,hh,Bfield,curlB,gradB=gradB)
        do i=1,npart
           B2i = dot_product(Bfield(:,i),Bfield(:,i))
           if (B2i > 1.e-8) then
              ! Tricco & Price (2013) resistivity switch
              alpha(3,i) = min(hh(i)*norm2(gradB(:,:,i))/sqrt(B2i + epsilon(B2i)),1.0)
           else
              alpha(3,i) = 0.
           endif
           !if (B2i > 1.e-8) then
              !call curl3D_epsijk(gradB(:,:,i),curlBi)
              !print*, ' gradB = ',gradB(:,:,i)
              !print*,' got curlz = ',curlBi(:)
              !print*,' should be = ',curlB(:,i)
              !read*
           !endif
        enddo
     else
        !--get J using standard (differenced) curl operator
        call get_curl(1,npart,x,pmass,rho,hh,Bfield,curlB)
     endif
  endif
!
!--calculate thermal energy from the conserved energy (or entropy)
!

  select case(iener)
  case(3)     ! total energy is evolved
     do i=1,npart
        v2i = DOT_PRODUCT(vel(:,i),vel(:,i))
        B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
        uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        if (uu(i).lt.0.) then
           nerr = nerr + 1
           uu(i) = 0.
        endif
     enddo
     if (nerr.gt.0) write(iprint,*) 'Warning: utherm -ve on ',nerr,' particles '
  case(1)  ! en = entropy variable
     uu = en/(gamma-1.)*rho**(gamma-1.)
  case(4)  ! en = rho*u (volume thermal energy variable)
     uu = en/rho
  case(5)
     call smooth_variable(en,uu,x,pmass,hh,rho)
  case default    ! en = thermal energy
     if (damp.gt.0.) then
        uu = en
     !   do i=1,npart
     !      dx = x(2,i) - xmin(2)
     !      uu(i) = przero/((gamma - 1.)*yprofile(dx,denszero,densmedium,smoothl))
     !      en(i) = uu(i)
     !   enddo
        !if (mod(nsteps,100).eq.0) then
        !   vel = 0.
        !   velin = 0.
        !endif
     else
        uu = en
     endif
  end select
!
!--get lambda for linear combination of kernels
!
 if (ibiascorrection.gt.0) then
    call get_lambda(sqrtg,ntotal)
 endif
!
!--get second derivatives of velocity
!
 if (.false.) then
    call get_2ndderivs(2,npart,x,pmass,rho,hh,vel,del2v,graddivv)
 endif
!
!--call equation of state calculation
!
!
!--compute unity (stored as psi) for later use in equation of motion
!  (use pressure temporarily here)
 if (iprterm.eq.12) then
    pr(:) = 1.0
    call smooth_variable(pr,psi,x,pmass,hh,rho)
    do i=1,npart
       if (psi(i).lt.0.999 .or. psi(i).gt.1.001) print*,'unity = ',i,x(1,i),psi(i)
    enddo
 endif

 if (iprterm.eq.10) then
    pr(1:npart) = (gamma-1.)*psi(1:npart)
    spsound(1:npart) = gamma*pr(1:npart)/rho(1:npart)
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) &
       stop 'iprterm=10 not implemented with one-fluid dust'
 elseif (iprterm.eq.11) then
    call equation_of_state(pr(1:npart),spsound(1:npart),uu(1:npart),  &
                        rho(1:npart),psi(1:npart))
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) &
       stop 'iprterm=11 not implemented with one-fluid dust'
 else
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       !
       !--for one fluid dust dens(i) is the GAS density, while rho(i) is the TOTAL density
       !
       call equation_of_state(pr(1:npart),spsound(1:npart),uu(1:npart),dens(1:npart))    
    else
       call equation_of_state(pr(1:npart),spsound(1:npart),uu(1:npart),rho(1:npart))
    endif
 endif
!
!--make fixed particles exact replicas of their closest particle
!
  if (any(ibound.eq.1)) then
     do i=1,npart
        if (itype(i).eq.itypebnd) then
           j = ireal(i)
           if (j > 0) call copy_particle(i,j)
           !vel(:,i) = vel(:,j)
           !where (ibound.eq.1) vel(:,i) = -vel(:,j)
!!           uu(i) = uu(j)
!           pr(i) = pr(j)
!           spsound(i) = spsound(j)
        endif
     enddo
  endif
!
!--copy the primitive variables onto the ghost particles
! 
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        curlB(:,i) = curlB(:,j)
        if (allocated(zeta)) zeta(i) = zeta(j)
        psi(i) = psi(j)
        if (allocated(dens)) dens(i) = dens(j)
        if (any(ibound.eq.6)) then
           pri = 2.5 - 0.1*dens(i)*x(2,i)
           uu(i) = pri/((gamma-1.)*dens(i))
           spsound(i) = sqrt(gamma*pri/dens(i))
        else
           uu(i) = uu(j)
           spsound(i) = spsound(j)
           pr(i) = pr(j)
        endif
        Bfield(:,i) = Bfield(:,j)
        if (idust > 0 .and. idust /= 2) then
           dustfrac(i) = dustfrac(j)
        endif
        if (all(ibound.eq.3)) call copy_particle(i,j) ! just to be sure
     enddo
  endif
  
  return
end subroutine conservative2primitive

!---------------------------------------------------------------------
! this subroutine is called after setting up the initial conditions
! to set the initial values of the conserved variables
! makes a call to density to calculate the density by summation
! also sets the value of the initial smoothing length
!---------------------------------------------------------------------
subroutine primitive2conservative
  use dimen_mhd
  use debug
  use loguns
  use bound, only:ireal
  use hterms, only:rhomin,gradgradh,gradh,zeta
  use eos
  use options
  use part
  use setup_params
  use timestep
  use getcurl
  use getBeulerpots, only:get_B_eulerpots
  use rates, only:gradpsi
  use derivB, only:curlB
  use utils,  only:minmaxave
!  use resistivity, only:Bdiffusion
  implicit none
  integer :: i,j,iktemp,iremap,nremap
  real :: B2i, v2i, hmin, hmax, hav, polyki, gam1, pri, dhdrhoi, emagold,emag
  real, dimension(ndimV) :: Binti
  logical :: isetpolyk,remap

  if (trace) write(iprint,*) ' Entering subroutine primitive2conservative'
!
!--calculate conserved density (and initial smoothing length)
!
  isetpolyk = .false.
  do i=1,npart
     if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
        ! Note: for most of the one-fluid dust setups, we set dens to mean the
        ! total density. This is only fixed once conservative2primitive has 
        ! been called. Hence do NOT divide by (1 - eps) below so that 
        ! smoothing lengths can be guessed correctly
        rho(i) = dens(i)  !/(1. - dustfrac(i)) ! rho is total mass density, dens is gas density only  
     else
        rho(i) = dens(i)
     endif
     if (usenumdens) then
        hh(i) = hfact*psep
     else
        hh(i) = hfact*(pmass(i)/(rho(i) + rhomin))**dndim
     endif
!
!--also work out what polyk should be if using iener = 0
!  (only do this once, otherwise give an error)
!
     if (iener.eq.0 .and. itype(i).eq.itypegas .and. iexternal_force.ne.10) then
        gam1 = gamma - 1.
        if (gamma.le.1.0001) then
           polyki = 2./3.*uu(i)
        else
           polyki = gam1*uu(i)/dens(i)**(gam1)
        endif
        if (abs(polyki-polyk)/polyk.gt.1.e-8) then
           write(iprint,*) 'NOTE: setting polyk = ',polyki,' (infile says ',polyk,')'
           if (isetpolyk) then
              write(iprint,*) 'ERROR: particle ',i,': multiple polyk factors in setup but using iener = 0'
              write(iprint,*) 'uu = ',uu(i),' dens = ',dens(i),' gamma = ',gamma,' polyk = ',gam1*uu(i)/dens(i)**(gam1),polyki
              stop
           endif
           isetpolyk = .true.
           polyk = polyki
        endif
     endif
  enddo
  if (ihvar.le.0) then
     call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
     hh(1:ntotal) = hav
     print*,'HMIN,hmax,hav = ',hmin,hmax,hav
  endif

!
!--overwrite this with a direct summation
!  
  if (icty.eq.0 .or. ndirect.lt.100000) then
     if (any(hh(1:npart).le.tiny(hh))) stop 'h < 0 in primitive2conservative'
     write(iprint,*) 'Calculating initial density...' 
     if (any(ibound > 1)) call set_ghost_particles
     call set_linklist
     iktemp = ikernav
     ikernav = 3                ! consistent with h for first density evaluation
     call iterate_density        ! evaluate density by direct summation
!     call densityiterate
     ikernav = iktemp  
!!     hh(1:ntotal) = hfact*(pmass(1:ntotal)/(rho(1:ntotal)+rhomin))**dndim
     if (ihvar.le.0) then
        call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
        hh(1:npart) = hav
     endif
  endif

  if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
     if (use_sqrtdustfrac) then
        dustevol = sqrt(dustfrac*rho)
     else
        dustevol = dustfrac
     endif
     dens = rho*(1. - dustfrac)
  else
     dens = rho
  endif
!
!--calculate conserved variable from the magnetic flux density B
!  
  remap = .false.
  x0 = x     ! setup initial positions for Jacobian mapping for all particles
  rho0 = rho ! initial densities for checking the Jacobian mapping

  select case(imhd)
  case(11:)    ! if using B or remapped B as conserved variable
     Bevol = Bfield
     if (imhd.eq.20) write(iprint,*) ' Using remapped B as conserved variable...'
  case(10)   ! remapped B/rho as conserved variable
     write(iprint,*) ' Using remapped B/rho as conserved variable...'
     do i=1,npart
        Bevol(:,i) = Bfield(:,i)/rho(i)
        !if (i.le.10) print*,i,' B/rho (before) = ',Bevol(:,i)
     enddo
     !--check that remapping at time zero does nothing
     !call get_B_eulerpots(3,npart,x,pmass,rho,hh,Bevol,x0,Bevol,remap=.true.)
     !do i=1,10
     !   print*,i,'B/rho (after) = ',Bevol(:,i)
     !enddo
     !read*
  case(1:9)   ! if using B/rho as conserved variable
     do i=1,npart
        Bevol(:,i) = Bfield(:,i)/rho(i)
     enddo
  case(-1,-2) ! if using vector potential
     if (iprterm.ge.10) stop 'conflict with use of psi variable (iprterm=10/imhd<0)'
     write(iprint,*) 'getting B field from vector potential (init)...'
     if (iuse_exact_derivs.gt.0) then
        !--get an exact linear curl A (iderivtype=2 input to get_B_eulerpots)
        call get_B_eulerpots(2,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap)
     else
        call get_curl(1,npart,x,pmass,rho,hh,Bevol,Bfield,gradpsi)
     endif
     do i=1,npart
        Binti(:) = Bfield(:,i)
        Bfield(:,i) = Bfield(:,i) + Bconst(:)
        !--construct xi term
        if (imagforce.eq.7 .and. ihvar.gt.0) then
           dhdrhoi = -hh(i)/(ndim*(rho(i)))
           gradpsi(:,i) = dhdrhoi*gradpsi(:,i)
           zeta(i) = dot_product(Bfield(:,i),gradpsi(:,i)) + dot_product(Bfield(:,i),Binti(:))*gradgradh(i)*gradh(i)
           psi(i) = dot_product(Bfield(:,i),Binti(:))*gradh(i)/rho(i)*dhdrhoi
           !print*,i,'xi = ',psi(i),zeta(i),' zeta = ',gradgradh(i),' gradpsi = ',gradpsi(:,i) !,Bfield(:,i)
        else
           zeta(i) = 0.
           psi(i) = 0.
        endif
     enddo
     !--reset gradpsi to zero after we have finished using it
     gradpsi(:,:) = 0.
     write(iprint,*) 'getting J from B (init)...'
     !--copy Bfield onto ghosts
     if (any(ibound.eq.1)) then
        do i=1,npart
           if (itype(i).eq.itypebnd) then
              j = ireal(i)
              if (j > 0) Bfield(:,i) = Bfield(:,j)
           endif
        enddo
     endif
     if (any(ibound.gt.1)) then
        do i=npart+1,ntotal
           j = ireal(i)
           if (j > 0) Bfield(:,i) = Bfield(:,j)
        enddo
     endif
     call get_curl(Jcurltype,npart,x,pmass,rho,hh,Bfield,curlB)
     !print*,' curlB(100) = ',curlB(:,100)
  case(-3) ! Generalised Euler potentials
     !--compute B = \nabla alpha_k x \nabla \beta^k
     write(iprint,*) 'getting B field from Generalised Euler Potentials (init)... '
     nremap = 1
     remap = (nremap.gt.1)
     
     do iremap=1,nremap
        call get_B_eulerpots(1,npart,x,pmass,rho,hh,Bevol,x0,Bfield,remap)

        if (remap) then
           !--remap x0 for all particles
           x0 = x
           !--copy remapped Bevol onto boundary/ghost particles
           if (any(ibound.eq.1)) then
              do i=1,npart
                 if (itype(i).eq.itypebnd) then ! fixed particles
                    j = ireal(i)
                    Bevol(:,i) = Bevol(:,j)
                 endif
              enddo
           endif
           if (any(ibound.gt.1)) then  ! ghost particles
              do i=npart+1,ntotal
                 j = ireal(i)
                 Bevol(:,i) = Bevol(:,j)
              enddo
           endif
        !
        !--recompute B field with remapped potentials
        !
           emag = emag_calc(pmass,rho,Bfield,npart)
           if (iremap.eq.1) then
              emagold = emag
              print*,' magnetic energy before remapping = ',emagold
           else
              print*,' magnetic energy after ',iremap,'nd/th remapping = ',emag, ' change = ',(emag-emagold)/emagold
           endif
        endif
     enddo
     
     !--copy Bfield onto ghosts
     if (any(ibound.eq.1)) then
        do i=1,npart
           if (itype(i).eq.itypebnd) then
              j = ireal(i)
              if (j > 0) Bfield(:,i) = Bfield(:,j)
           endif
        enddo
     endif
     if (any(ibound.gt.1)) then
        do i=npart+1,ntotal
           j = ireal(i)
           Bfield(:,i) = Bfield(:,j)
        enddo
     endif

  end select
!
!--calculate conserved energy (or entropy) from the thermal energy
!
  select case(iener)
  case(3)     ! total energy is evolved
     do i=1,npart
        v2i = DOT_PRODUCT(vel(:,i),vel(:,i))
        B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
        en(i) = uu(i) + 0.5*v2i + 0.5*B2i
        if (uu(i).lt.0.) stop 'primitive2conservative: utherm -ve '
     enddo
  case(1)   ! en = entropy
     where (dens > 0.)
        en = (gamma-1.)*uu/dens**(gamma-1.)
     end where
  case(4)   ! en = rho*u (volume thermal energy)
     en = uu*dens
  case default        ! en = thermal energy
     en = uu
  end select
!
!--call equation of state calculation
!  (not ghosts, but including fixed particles)
  if (iprterm.eq.11) then
     call equation_of_state(pr(1:npart),spsound(1:npart), &
                            uu(1:npart),dens(1:npart),psi(1:npart))
  else
     call equation_of_state(pr(1:npart),spsound(1:npart), &
                            uu(1:npart),dens(1:npart))
  endif
!
!--copy the conservative variables onto the ghost particles
!  
  if (any(ibound > 1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        if (j > 0) then
           call copy_particle(i,j)
           if (any(ibound.eq.6)) then
              pri = 2.5 - 0.1*dens(i)*x(2,i)
              uu(i) = pri/((gamma-1.)*dens(i))
              spsound(i) = sqrt(gamma*pri/dens(i))
              en(i) = uu(i)
              if (iener.ne.2) stop 'iener.ne.2 not implemented for ibound=6'
           else
              en(i) = en(j)
              pr(i) = pr(j)
              spsound(i) = spsound(j)
           endif
        endif
     enddo
  endif
!
!--make fixed particles exact replicas of their closest particle
!
  if (any(ibound.eq.1) .and. imhd.lt.0) then
     do i=1,npart
        if (itype(i).eq.itypebnd) then
           j = ireal(i)
           if (j > 0) then
              call copy_particle(i,j)
              vel(:,i) = vel(:,j)
           endif
        endif
     enddo
  endif

!
!--call rates to get initial timesteps, div B etc
!
  call derivs
 
  return  
end subroutine primitive2conservative

real function emag_calc(pmass,rho,Bfield,npart)
  implicit none
  integer, intent(in) :: npart
  real, dimension(:), intent(in) :: pmass,rho
  real, dimension(:,:), intent(in) :: Bfield
  real :: emag, B2i
  integer :: i
  
  emag = 0.
  do i=1,npart
     B2i = dot_product(Bfield(:,i),Bfield(:,i))
     emag = emag + pmass(i)*0.5*B2i/rho(i)
  enddo
  emag_calc = emag

end function emag_calc

end module cons2prim
