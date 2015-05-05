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
  use rates,  only:gradpsi
  use hterms, only:gradgradh,gradh,zeta
  use derivB, only:curlB
  use part_in, only:dustevolin !,velin
  implicit none
  integer :: i,j,nerr,nerr1
  real :: B2i, v2i, pri, dhdrhoi
  real, dimension(ndimV) :: Binti,Bfieldi
  logical, parameter :: JincludesBext = .true.

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
  if (nerr > 0) print*,'WARNING: dust fraction < 0 on ',nerr,' particles'
  if (nerr1 > 0) print*,'WARNING: dust fraction > 1 on ',nerr1,' particles'
!
!--calculate magnetic flux density B from the conserved variable
!  
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bfield = Bevol
  elseif (imhd.gt.0) then ! if using B/rho as conserved variable
     do i=1,npart
        Bfield(:,i) = Bevol(:,i)*rho(i)
     enddo
  elseif (imhd.lt.0) then ! if using vector potential
     !write(iprint,*) 'getting B field from vector potential...'
     if (iprterm.ge.10) stop 'conflict with use of psi variable (iprterm=10/imhd<0)'
     gradpsi(:,:) = 0.
     psi(:) = 0.
     zeta(:) = 0.
     call get_curl(1,npart,x,pmass,rho,hh,Bevol,Bfield,gradpsi)
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
  endif
!
!--get magnetic current, needed for ambipolar diffusion calculation
!
  if (iambipolar > 0) then
     !--get J using standard (differenced) curl operator
     call get_curl(1,npart,x,pmass,rho,hh,Bfield,curlB)
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
  case default    ! en = thermal energy
     if (damp.gt.0.) then
        uu = en
        !if (mod(nsteps,100).eq.0) then
        !   vel = 0.
        !   velin = 0.
        !endif
     else
        uu = en
     endif
  end select
!
!--call equation of state calculation
!
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
  use rates, only:gradpsi
  use derivB, only:curlB
  use utils,  only:minmaxave
  implicit none
  integer :: i,j,iktemp
  real :: B2i, v2i, hmin, hmax, hav, polyki, gam1, pri, dhdrhoi
  real, dimension(ndimV) :: Binti
  logical :: isetpolyk

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
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bevol = Bfield
  elseif (imhd.gt.0) then ! if using B/rho as conserved variable
     do i=1,npart
        Bevol(:,i) = Bfield(:,i)/rho(i)
     enddo
  elseif (imhd.lt.0) then ! if using vector potential
     if (iprterm.ge.10) stop 'conflict with use of psi variable (iprterm=10/imhd<0)'
     write(iprint,*) 'getting B field from vector potential (init)...'
     call get_curl(1,npart,x,pmass,rho,hh,Bevol,Bfield,gradpsi)
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
  endif
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

end module cons2prim
