!---------------------------------------------------------------------
! these subroutines calculates the primitive variables from the
! conservative variables and vice versa.
!
! This version is for non-relativistic, cylindrical hydro (no MHD yet).
! The thermal energy is calculated from the
! total energy, whilst the magnetic field B is calculated from
! the conserved variable B/rho.
!
! also returns the source terms for the momentum equations
! and the determinant of the metric sqrtg
!
! These subroutines would be used in the GR case, for all the 
! variables. The evaluation would be far more complicated, however.
!---------------------------------------------------------------------
subroutine conservative2primitive(rho,vel,uu,en,Bfield,Bcons,ierr)
  use dimen_mhd
  use options
  use loguns  
  implicit none
  integer, intent(out) :: ierr  ! error code
  real, intent(in) :: en,rho
  real, intent(in), dimension(ndimV) :: vel,Bcons 
  real, intent(out) :: uu
  real, intent(out), dimension(ndimV) :: Bfield
  real :: B2i, v2i

  ierr = 0
!
!--calculate co-ordinate velocities (x derivatives) from the
!  conserved momentum variable
!
  select case (igeom)
  case(1)   ! cylindrical
     sqrtg = x(1)
     dens = rho/sqrtg
     vel(:) = pmom(:)
     if (ndimV.ge.2) vel(2) = pmom(2)/x(1)
  case(2)   ! spherical
     if (ndim.eq.1) then
        sqrtg = x(1)**2
     else
        sqrtg = (x(1)**2)*SIN(x(2))
     endif
     vel(1) = pmom(1)
     vel(2) = pmom(2)/x(1)
     vel(3) = pmom(3)/(x(1)*
  case (default) ! cartesian
     sqrtg = 1.
     dens = rho
     vel(:) = pmom(:)
  if (iener.eq.3) then     ! total energy is evolved
     v2i = DOT_PRODUCT(vel,vel)
     uu = en - 0.5*v2i - 0.5*B2i
  elseif (iener.ge.1) then ! thermal energy is evolved
     uu = en
  endif

  end select

!
!--calculate magnetic flux density B from the conserved variable
!  
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bfield(:) = Bcons(:)
     B2i = DOT_PRODUCT(Bcons,Bcons)/rho
  elseif (imhd.ne.0) then ! if using B/rho as conserved variable
     Bfield(:) = Bcons(:)*rho
     B2i = DOT_PRODUCT(Bcons,Bcons)*rho
  else
     B2i = 0.
  endif

  if (uu.lt.0.) then
     write(iprint,*) 'Warning: utherm -ve'
     ierr = 1  ! error if negative thermal energy
  endif

  return
end subroutine conservative2primitive

!---------------------------------------------------------------------
! this subroutine is called after setting up the initial conditions
! to set the initial values of the conserved variables
!---------------------------------------------------------------------
subroutine primitive2conservative(rho,vel,uu,en,Bfield,Bcons,ierr)
  use dimen_mhd
  use options
  implicit none
  integer, intent(out) :: ierr ! error code
  real, intent(out) :: en
  real, intent(out), dimension(ndimV) :: Bcons 
  real, intent(in) :: uu,rho
  real, intent(in), dimension(ndimV) :: vel,Bfield
  real :: B2i, v2i

  ierr = 0
!
!--calculate conserved momentum variable from 
!  coordinate velocities (x derivatives)
!
  select case(igeom)
  case(1)     ! cylindrical
     pmom(:) = vel(:)
     if (ndimV.ge.2) pmom(2) = x(1)*vel(2)
  case default
     pmom(:) = vel(:)
  end select
!
!--calculate conserved variable from the magnetic flux density B
!  
  if (imhd.ge.11) then    ! if using B as conserved variable
     Bcons(:) = Bfield(:)
     B2i = DOT_PRODUCT(Bcons,Bcons)/rho
  elseif (imhd.ne.0) then ! if using B/rho as conserved variable
     Bcons(:) = Bfield(:)/rho
     B2i = DOT_PRODUCT(Bcons,Bcons)*rho
  else
     B2i = 0.
  endif
!
!--calculate conserved energy (or entropy) from the thermal energy
!
  if (iener.eq.3) then     ! total energy is evolved
     v2i = DOT_PRODUCT(vel,vel)
     en = uu + 0.5*v2i + 0.5*B2i
     if (uu.lt.0.) ierr = 1  ! error if negative thermal energy
  elseif (iener.ge.1) then ! thermal energy is evolved
     en = uu
  endif

  return  
end subroutine primitive2conservative
