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
  implicit none
  integer :: i,j
  real :: B2i, v2i

  if (trace) write(iprint,*) ' Entering subroutine conservative2primitive'

  sqrtg = 1.
  dens = rho
!
!--calculate magnetic flux density B from the conserved variable
!  
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bfield = Bcons
  elseif (imhd.ne.0) then ! if using B/rho as conserved variable
     do i=1,npart
        Bfield(:,i) = Bcons(:,i)*rho(i)
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
           !write(iprint,*) 'Warning: utherm -ve, particle ',i
           uu(i) = 0.
        endif
     enddo
  else		 ! en = thermal energy
     uu = en
  endif
!
!--call equation of state calculation
!
  do i=1,npart
     call equation_of_state(pr(i),spsound(i),uu(i),rho(i),gamma)
  enddo
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
  use loguns
  use dimen_mhd
  use eos
  use options
  use part
  use setup_params
  use timestep
  implicit none
  integer :: i,iktemp
  real :: B2i, v2i
!
!--calculate conserved density (and initial smoothing length
!
  rho = dens
  hh(1:npart) = hfact*(pmass(1:npart)/rho(1:npart))**hpower
!
!--overwrite this with a direct summation
!  
  if (icty.eq.0 .or. ndirect.lt.100000) then
     write(iprint,*) 'Calculating initial density...' 
     if (ANY(ibound.GT.1)) call set_ghost_particles
     call link
     iktemp = ikernav
     ikernav = 3		! consistent with h for first density evaluation
     call iterate_density	! evaluate density by direct summation
     ikernav = iktemp  
     hh(1:npart) = hfact*(pmass(1:npart)/rho(1:npart))**hpower
  endif
!
!--calculate conserved variable from the magnetic flux density B
!  
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bcons = Bfield
  elseif (imhd.ne.0) then ! if using B/rho as conserved variable
     do i=1,npart
        Bcons(:,i) = Bfield(:,i)/rho(i)
     enddo
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
  else		! en = thermal energy
     en = uu
  endif
!
!--call equation of state calculation
!
  do i=1,npart	! not ghosts, but including fixed particles
    call equation_of_state(pr(i),spsound(i),uu(i),dens(i),gamma)     
  enddo
  
  return  
end subroutine primitive2conservative
