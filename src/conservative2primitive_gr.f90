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
  use eos
  use part

  use grutils
  implicit none
  real, dimension(ndimV) :: gdiag
  real :: B2i, v2i
!
!--calculate co-ordinate velocities (x derivatives) from the
!  conserved momentum variable
!
  sqrtg = 1.

  select case (igeom)
!
!--cylindrical metric
!
  case(1)
     do i=1,npart
        !--metric terms
        gdiag(:) = 1.
        if (ndimV.ge.2) gdiag(2) = x(1,i)**2
        sqrtg(i) = abs(x(1,i))

        dens(i) = rho(i)/sqrtg(i)
        vel(:,i) = pmom(:,i)
        if (ndimV.ge.2) vel(2,i) = pmom(2,i)/(x(1,i)**2)

        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
           B2i = 0. ! no mhd yet
           uu(i) = en(i) - 0.5*v2i - 0.5*B2i
           if (uu(i).lt.0.) then
              write(iprint,*) 'Warning: utherm -ve, particle ',i
              uu(i) = 0.
           endif
        else   ! thermal energy is evolved
           uu(i) = en(i)
        endif
        !
        !--call equation of state to get pressure
        !
        call equation_of_state(pr(i),vsound(i),uu(i),dens(i),gamma)
        !!--compute source terms (derivatives of metric) for momentum equation
        sourceterms(1,i) = pr(i)/rho(i)
        if (ndimV.ge.2) sourceterms(1,i) = sourceterms(1,i) + pmom(2,i)
     enddo
!
!--spherical polars
!
  case(2)
     do i=1,npart
        !--metric terms
        gdiag(:) = 1.
        if (ndimV.ge.2) gdiag(2) = x(1,i)**2
        if (ndimV.ge.3) gdiag(3) = gdiag(2)*SIN(x(2,i))
        if (ndim.eq.1) then
           sqrtg(i) = x(1,i)**2
        else
           sqrtg(i) = SIN(x(2,i))*x(1,i)**2
        endif

        dens(i) = rho(i)/sqrtg(i)
        vel(1,i) = pmom(1,i)
        if (ndimV.ge.2) vel(2,i) = pmom(2,i)/(x(1,i)**2)
        if (ndimV.eq.3) vel(3,i) = pmom(3,i)/sqrtg(i)
        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
           uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        else       ! thermal energy is evolved
           uu(i) = en(i)
        endif
        !
        !--call equation of state to get pressure
        !
        call equation_of_state(pr(i),vsound(i),uu(i),dens(i),gamma)
        !!--compute source terms (derivatives of metric) for momentum equation
        sourceterms(1,i) = pr(i)/rho(i)
        if (ndimV.ge.2) sourceterms(1,i) = sourceterms(1,i) + x(1,i)**2*vel(2,i)
        if (ndimV.ge.3) sourceterms(2,i) = x(1,i)**2*vel(3,i)
     enddo
!
!--spherical polars with logarithmic radial co-ordinate
!
  case(3)
     do i=1,npart
        !--metric terms
        gdiag = 1.
        
        rr = exp(x(1,i))
        if (ndim.eq.1) then
           sqrtg(i) = rr**3
        else
           sqrtg(i) = SIN(x(2,i))*rr**3
        endif
        vel(1,i) = pmom(1,i)
        if (ndimV.ge.2) vel(2,i) = pmom(2,i)/(x(1,i)**2)
        if (ndimV.ge.3) vel(3,i) = pmom(3,i)/sqrtg(i)
        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT(vel,vel)
           uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        else       ! thermal energy is evolved
           uu(i) = en(i)
        endif
        sourceterms(:,i) = 0.!! WORK THESE OUT!!
     enddo
!
!--Minkowski metric (special relativity)
!
  case(10)
     stop 'special rel not yet implemented'
!
!--cartesian co-ordinates
!
  case (default)
     dens = rho
     !!vel = pmom   
     !
     !--calculate magnetic flux density B from the conserved variable
     !  
     if (imhd.ge.11) then    ! if using B as conserved variable
        Bfield = Bcons
     elseif (imhd.ne.0) then ! if using B/rho as conserved variable
        Bfield = Bcons*rho
     endif
     do i=1,npart
        if (iener.eq.3) then      ! total energy is evolved
           v2i = DOT_PRODUCT(vel,vel)
           B2i = DOT_PRODUCT(Bfield,Bfield)/rho(i)
           uu(i) = en(i) - 0.5*v2i - 0.5*B2i
           if (uu(i).lt.0.) then
              write(iprint,*) 'Warning: utherm -ve, particle ',i
              uu(i) = 0.
           endif
        else
           uu(i) = en(i)
        endif
        !
        !--call equation of state to get pressure
        !
        call equation_of_state(pr(i),spsound(i),uu(i),dens(i),gamma)           
     enddo
           
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
  use options
  implicit none
  real, dimension(ndim) :: gdiag
  real :: B2i, v2i
!
!--calculate conserved momentum variable from 
!  coordinate velocities (x derivatives)
!
  sqrtg = 1.

  select case(igeom)
!
!--cylindrical co-ords
!
  case(1)
     do i=1,npart
        !--metric terms
        gdiag(:) = 1.
        if (ndimV.ge.2) gdiag(2) = x(1,i)**2
        sqrtg(i) = abs(x(1,i))
        rho(i) = dens(i)*sqrtg(i)
        pmom(:,i) = vel(:,i)
        if (ndimV.ge.2) pmom(2,:) = x(1,:)*vel(2,:)
        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT_GR(vel(:,i),vel(:,i),gdiag,ndimV)
           B2i = 0. ! no mhd yet
           en(i) = uu(i) + 0.5*v2i + 0.5*B2i
           if (uu(i).lt.0.) then
              write(iprint,*) 'Warning: utherm -ve, particle ',i
              uu(i) = 0.
           endif
        else   ! thermal energy is evolved
           uu(i) = en(i)
        endif
     enddo
!
!--spherical co-ords
!
  case(2)
     if (ndim.eq.1) then
        sqrtg(:) = x(1,:)**2
     else
        sqrtg(:) = SIN(x(2,:))*x(1,:)**2
     endif
     rho = dens*sqrtg
     do i=1,npart
        vel(1,i) = pmom(1,i)
        if (ndimV.ge.2) vel(2,i) = pmom(2,i)/(x(1,i)**2)
        if (ndimV.eq.3) vel(3,i) = pmom(3,i)/sqrtg(i)
        if (iener.eq.3) then     ! total energy is evolved
           v2i = DOT_PRODUCT(vel(:,i),vel(:,i))
           uu(i) = en(i) - 0.5*v2i - 0.5*B2i
        else       ! thermal energy is evolved
           uu(i) = en(i)
        endif
     enddo
!
!--cartesian
!
  case default
     pmom(:) = vel(:)
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
        en = uu(i) + 0.5*v2i + 0.5*B2i
        if (uu(i).lt.0.) ierr = 1  ! error if negative thermal energy
     elseif (iener.ge.1) then ! thermal energy is evolved
        en = uu(i)
     endif
  end select

  return  
end subroutine primitive2conservative
