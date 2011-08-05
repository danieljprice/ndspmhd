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
subroutine conservative2primitive
  use dimen_mhd
  use debug
  use options
  use loguns
  use bound
  use eos
  use part
  use getcurl
  implicit none
  integer :: i,j,nerr
  real :: B2i, v2i

  if (trace) write(iprint,*) ' Entering subroutine conservative2primitive'

  nerr = 0
  sqrtg = 1.
  dens = rho
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
     write(iprint,*) 'getting B field from vector potential...'
     call get_curl(npart,x,pmass,rho,hh,Bevol,Bfield)
     do i=1,npart
        Bfield(:,i) = Bfield(:,i) + Bconst(:)
     enddo
  endif
!
!--calculate thermal energy from the conserved energy (or entropy)
!
  if (iener.eq.3) then     ! total energy is evolved
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
  elseif (iener.eq.1) then  ! en = entropy variable
     uu = en/(gamma-1.)*rho**(gamma-1.)
  else                 ! en = thermal energy
     uu = en
  endif
!
!--call equation of state calculation
!
 call equation_of_state(pr(1:npart),spsound(1:npart),uu(1:npart),  &
                        rho(1:npart),gamma,polyk,npart)
!
!--make fixed particles exact replicas of their closest particle
!
  if (any(ibound.eq.1)) then
     do i=1,npart
        if (itype(i).eq.1) then
           j = ireal(i)
           uu(i) = uu(j)
           pr(i) = pr(j)
           spsound(i) = spsound(j)
           Bfield(:,i) = Bfield(:,j)
        endif
     enddo
  endif
!
!--copy the primitive variables onto the ghost particles
! 
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        if (allocated(dens)) dens(i) = dens(j)
        uu(i) = uu(j)
        spsound(i) = spsound(j)
        pr(i) = pr(j)
        Bfield(:,i) = Bfield(:,j)
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
  use hterms, only:rhomin
  use eos
  use options
  use part
  use setup_params
  use timestep
  use getcurl
  implicit none
  integer :: i,j,iktemp
  real :: B2i, v2i, hmin, hmax, hav, polyki, gam1
  logical :: isetpolyk

  if (trace) write(iprint,*) ' Entering subroutine primitive2conservative'
!
!--calculate conserved density (and initial smoothing length)
!
  isetpolyk = .false.
  do i=1,npart
     rho(i) = dens(i)
     hh(i) = hfact*(pmass(i)/(rho(i) + rhomin))**dndim
!
!--also work out what polyk should be if using iener = 0
!  (only do this once, otherwise give an error)
!
     if (iener.eq.0) then
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
  endif
!
!--overwrite this with a direct summation
!  
  if (icty.eq.0 .or. ndirect.lt.100000) then
     write(iprint,*) 'Calculating initial density...' 
     if (ANY(ibound.GT.1)) call set_ghost_particles
     call set_linklist
     iktemp = ikernav
!     ikernav = 3                ! consistent with h for first density evaluation
     call iterate_density        ! evaluate density by direct summation
     ikernav = iktemp  
!!     hh(1:ntotal) = hfact*(pmass(1:ntotal)/(rho(1:ntotal)+rhomin))**dndim
     if (ihvar.le.0) then
        call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
        hh(1:npart) = hav
     endif
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
     write(iprint,*) 'getting B field from vector potential (init)...'
     call get_curl(npart,x,pmass,rho,hh,Bevol,Bfield)
  endif
!
!--calculate conserved energy (or entropy) from the thermal energy
!
  if (iener.eq.3) then     ! total energy is evolved
     do i=1,npart
        v2i = DOT_PRODUCT(vel(:,i),vel(:,i))
        B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
        en(i) = uu(i) + 0.5*v2i + 0.5*B2i
        if (uu(i).lt.0.) stop 'primitive2conservative: utherm -ve '
     enddo
  elseif (iener.eq.1) then ! en = entropy
     en = (gamma-1.)*uu/rho**(gamma-1.)
  else                ! en = thermal energy
     en = uu
  endif
!
!--call equation of state calculation
!  (not ghosts, but including fixed particles)
  call equation_of_state(pr(1:npart),spsound(1:npart), &
                         uu(1:npart),dens(1:npart),gamma,polyk,npart)  
!
!--copy the conservative variables onto the ghost particles
!  
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        en(i) = en(j)
        pr(i) = pr(j)
        spsound(i) = spsound(j)
        Bevol(:,i) = Bevol(:,j)
        Bfield(:,i) = Bfield(:,j)
!        call copy_particle(i,j)
     enddo
  endif
!
!--call rates to get initial timesteps, div B etc
!
  call get_rates

  return  
end subroutine primitive2conservative
