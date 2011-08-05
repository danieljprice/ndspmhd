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
!
!--calculate thermal energy from the conserved energy (or entropy)
!
  if (iener.eq.3) then     ! total energy is evolved
     v2i = DOT_PRODUCT(vel,vel)
     uu = en - 0.5*v2i - 0.5*B2i
  else		 ! en = thermal energy
     uu = en
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
     if (uu.lt.0.) stop 'primitive2conservative: utherm -ve'
  else		! en = thermal energy
     en = uu
  endif

  return  
end subroutine primitive2conservative
