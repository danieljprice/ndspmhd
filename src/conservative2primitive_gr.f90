!---------------------------------------------------------------------
! these subroutines calculates the primitive variables from the
! conservative variables and vice versa.
!
! This version is for non-relativistic hydro in arbitrary co-ordinate
! systems (no MHD yet except in cartesian).
! The thermal energy is calculated from the
! total energy
! In cartesian coords the magnetic field B is calculated from
! the conserved variable B/rho.
!
! The equation of state is called to calculate the pressure
!
! also returns the source terms for the momentum equations
! and the determinant of the metric sqrtg
!
! This subroutine is more complicated in the relativistic case since the 
! Lorentz factor is involved everywhere.
!---------------------------------------------------------------------
subroutine conservative2primitive
  use dimen_mhd
  use options
  use loguns  
  use bound
  use eos
  use part
  use grutils
  implicit none
  integer :: i,j
  real, dimension(ndimV) :: gdiag
  real :: B2i, v2i

  sqrtg = 1.
!
!--this is the procedure for general (non-relativistic) diagonal metrics
!
  do i=1,npart
     call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,igeom)
     
     dens(i) = rho(i)/sqrtg(i)
     if (allocated(pmom)) vel(:,i) = pmom(:,i)/gdiag(:)
     
     if (imhd.ge.11) then    ! if using B as conserved variable
        Bfield(:,i) = Bcons(:,i)
     elseif (imhd.ne.0) then ! if using B/rho as conserved variable
        Bfield(:,i) = Bcons(:,i)*rho(i)
     endif     
     
     if (iener.eq.3) then     ! total energy is evolved
        v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
        B2i = DOT_PRODUCT_GR(Bfield(:,i),Bfield(:,i),gdiag,ndimV)/rho(i)
        uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        if (uu(i).lt.0.) then
           write(iprint,*) 'Warning: utherm -ve, particle ',i
           uu(i) = 0.
        endif
     else   ! thermal energy is evolved
        uu(i) = en(i)
     endif
     !
     !--call equation of state to get pressure (needed for source terms)
     !
     call equation_of_state(pr(i),spsound(i),uu(i),dens(i),gamma)
  enddo

!
!--compute source terms (derivatives of metric) for momentum equation
!  (these are 1/(2*dens)* T^\mu\nu d(g_\mu_\nu)/dx^i
!
  select case(igeom)
  !
  !--cylindrical
  !
  case(1)
     sourceterms(1,:) = pr(:)/rho(:)
     if (ndimV.ge.2) sourceterms(1,:) = sourceterms(1,:) + pmom(2,:)
  !
  !--spherical
  !
  case(2)
     do i=1,npart
        sourceterms(1,i) = pr(i)/rho(i)
        if (ndimV.ge.2) sourceterms(1,i) = sourceterms(1,i) + x(1,i)**2*vel(2,i)
        if (ndimV.ge.3) sourceterms(2,i) = x(1,i)**2*vel(3,i)
     enddo
  !
  !--spherical polars with logarithmic radial co-ordinate
  !
  case(3)
     do i=1,npart
        sourceterms(:,i) = 0.!! WORK THESE OUT!!
     enddo
  !
  !--note that source terms are zero and not allocated for cartesian
  !         
  end select
  
!
!--copy the primitive variables onto the ghost particles
! 
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        dens(i) = dens(j)
        vel(:,i) = vel(:,j)
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
!---------------------------------------------------------------------
subroutine primitive2conservative
  use dimen_mhd
  use eos
  use options
  use part
  use grutils
  use setup_params
  use timestep
  implicit none
  integer :: i,iktemp
  real, dimension(ndim) :: gdiag
  real :: B2i, v2i

  sqrtg = 1.
!
!--this is the procedure for general (non-relativistic) diagonal metrics
!
  if (igeom.lt.10) then
     do i=1,npart
        call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,igeom)
        
        rho(i) = dens(i)*sqrtg(i)
	hh(i) = hfact*(pmass(i)/rho(i))**hpower
     enddo
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
	
     do i=1,npart
        call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,igeom)
	
        if (allocated(pmom)) pmom(:,i) = vel(:,i)*gdiag(:)
        
        if (imhd.ge.11) then    ! if using B as conserved variable
           Bcons(:,i) = Bfield(:,i)
        elseif (imhd.ne.0) then ! if using B/rho as conserved variable
           Bcons(:,i) = Bfield(:,i)/rho(i)
        endif
        
        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
           B2i = DOT_PRODUCT_GR(Bfield(:,i),Bfield(:,i),gdiag,ndimV)/rho(i)
           en(i) = uu(i) + 0.5*v2i + 0.5*B2i
        else   ! thermal energy is evolved
           en(i) = uu(i)
        endif
	
        call equation_of_state(pr(i),spsound(i),uu(i),rho(i)/sqrtg(i),gamma) 	
     enddo
  endif

  return  
end subroutine primitive2conservative
