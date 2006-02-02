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
  use debug
  use options
  use loguns  
  use bound
  use eos
  use part
  use grutils
  implicit none
  integer :: i,j,nerr
  real, dimension(ndimV) :: gdiag
  real :: B2i, v2i, dsqrtgi

  if (trace) write(iprint,*) ' Entering subroutine conservative2primitive'

  nerr = 0
  sqrtg = 1.
!
!--this is the procedure for general (non-relativistic) diagonal metrics
!
  do i=1,npart
     call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,geom)
     
     dsqrtgi = 1./sqrtg(i)
     dens(i) = rho(i)*dsqrtgi
     if (allocated(pmom)) vel(:,i) = pmom(:,i)/gdiag(:)
     
     if (imhd.ge.11) then    ! if using B* as conserved variable
        Bfield(:,i) = Bevol(:,i)*dsqrtgi
     elseif (imhd.gt.0) then ! if using B*/rho* as conserved variable
        Bfield(:,i) = Bevol(:,i)*dens(i)  ! this is rho(i)*sqrtg(i)
     elseif (imhd.lt.0) then
        stop 'GR vector potential not implemented'
     endif     
     
     if (iener.eq.3) then     ! total energy is evolved
        v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
        B2i = DOT_PRODUCT_GR(Bfield(:,i),Bfield(:,i),gdiag,ndimV)/dens(i)
        uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        if (uu(i).lt.0.) then
           nerr = nerr + 1
           uu(i) = 0.
        endif
        if (nerr.gt.0) write(iprint,*) 'Warning: utherm -ve on ',nerr,' particles '
     elseif (iener.eq.1) then ! en = entropy variable
        stop 'GR entropy not implemented'
     else   ! thermal energy is evolved
        uu(i) = en(i)
     endif
     !
     !--call equation of state to get pressure (needed for source terms)
     !
     call equation_of_state(pr(i),spsound(i),uu(i),dens(i),gamma,polyk,1)
  enddo

  if (geom(1:4).ne.'cart') then
     if (allocated(sourceterms)) then
        call metric_derivs(x(:,1:npart),vel(:,1:npart),rho(1:npart),pr(1:npart), &
                           sourceterms(:,1:npart),ndim,ndimV,npart,geom)
     else
        write(iprint,*) 'ERROR: non-cartesian geometry but source terms not allocated'
        stop
     endif
  endif
!
!--copy the primitive variables onto the ghost particles
! 
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        call copy_particle(i,j)
        where(ibound.eq.2)
           pmom(:,i) = -pmom(:,j)
           vel(:,i) = -vel(:,j)
        end where
     enddo
  endif
 
  return
end subroutine conservative2primitive

!---------------------------------------------------------------------
! this subroutine is a cut-down version of conservative2primitive
! which just calculates v from pmom (used in timestepping)
!---------------------------------------------------------------------
subroutine getv_from_pmom(xi,pmomi,veli)
 use dimen_mhd
 use options, only:geom
 use grutils, only:metric_diag
 real, dimension(ndim), intent(in) :: xi
 real, dimension(ndimV), intent(in) :: pmomi
 real, dimension(ndimV), intent(out) :: veli
 real, dimension(ndimV) :: gdiag
 real :: sqrtgi

 call metric_diag(xi,gdiag,sqrtgi,ndim,ndimV,geom)
 veli(:) = pmomi(:)/gdiag(:)

end subroutine getv_from_pmom

!---------------------------------------------------------------------
! this subroutine is called after setting up the initial conditions
! to set the initial values of the conserved variables
!---------------------------------------------------------------------
subroutine primitive2conservative
  use dimen_mhd
  use debug
  use loguns
  
  use bound, only:ireal
  use eos
  use options
  use part
  use grutils
  use setup_params
  use timestep
  implicit none
  integer :: i,j,iktemp
  real, dimension(ndimV) :: gdiag
  real :: B2i, v2i, hmin, hmax, hav

  if (trace) write(iprint,*) ' Entering subroutine primitive2conservative'

  sqrtg = 1.
!
!--this is the procedure for general (non-relativistic) diagonal metrics
!
  do i=1,npart
     call metric_diag(x(1:ndim,i),gdiag,sqrtg(i),ndim,ndimV,geom)
     rho(i) = dens(i)*sqrtg(i)
     hh(i) = hfact*(pmass(i)/rho(i))**dndim
  enddo
!
!--overwrite this with a direct summation
!  
  if (icty.eq.0 .or. ndirect.lt.100000) then
     write(iprint,*) 'Calculating initial density...' 
     if (ANY(ibound.GT.1)) call set_ghost_particles
     call set_linklist
     iktemp = ikernav
!        ikernav = 3		! consistent with h for first density evaluation
     call iterate_density	! evaluate density by direct summation
     ikernav = iktemp  
!!        hh(1:npart) = hfact*(pmass(1:npart)/rho(1:npart))**dndim
     if (ihvar.le.0) then
        call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
        hh(1:npart) = hav
     endif
  endif

  do i=1,npart
     call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,geom)

     if (allocated(pmom)) pmom(:,i) = vel(:,i)*gdiag(:)

     if (imhd.ge.11) then    ! if using B* as conserved variable
        Bevol(:,i) = Bfield(:,i)*sqrtg(i)
     elseif (imhd.gt.0) then ! if using B*/rho* as conserved variable
        Bevol(:,i) = Bfield(:,i)*sqrtg(i)/rho(i)
     elseif (imhd.lt.0) then ! vector potential
        stop 'vector potential not implemented for GR'
     endif

     if (iener.eq.3) then     ! total energy is evolved
        v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
        B2i = DOT_PRODUCT_GR(Bfield(:,i),Bfield(:,i),gdiag,ndimV)/dens(i)
        en(i) = uu(i) + 0.5*v2i + 0.5*B2i
     else   ! thermal energy is evolved
        en(i) = uu(i)
     endif

     call equation_of_state(pr(i),spsound(i),uu(i),rho(i)/sqrtg(i), &
	  gamma,polyk,1) 	
  enddo
     
  if (geom(1:4).ne.'cart') then
     if (allocated(sourceterms)) then
        call metric_derivs(x(:,1:npart),vel(:,1:npart),rho(1:npart),pr(1:npart), &
                           sourceterms(:,1:npart),ndim,ndimV,npart,geom)
     else
        write(iprint,*) 'ERROR: non-cartesian geometry but source terms not allocated'
        stop
     endif
  endif
!
!--copy the conservative variables onto the ghost particles
!  
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        call copy_particle(i,j)
        where(ibound.eq.2)
           pmom(:,i) = -pmom(:,j)
           vel(:,i) = -vel(:,j)
        end where
     enddo
  endif
!
!--call rates to get initial timesteps, div B etc
!
  call get_rates     

  return  
end subroutine primitive2conservative
