!----------------------------------------------------------------
!     Set up for the 2D cartesian MRI test problem described in 
!     Hawley, Gammie & Balbus, (1995)
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i,iseed
 real :: massp,totmass
 real :: denszero,uuzero,cs0,polyk0
 real :: przero,Bzeroz,asize,betamhd,ran1,wavekmin,valfvenz
 real :: ksi,dksidX0,wavenum,wavenumsq,U,Omega2_wave
 integer, parameter :: ipert = 2
 integer, parameter :: isheartest = 1
 real, parameter    :: ksi0 = 1.d-5
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup (mri)'
 write(iprint,*) '2D magneto-rotational instability'
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 if (ndimV.ne.3) stop 'error: we need ndimV=3 for this problem'
!
!--geometry is cylindrical r-z
!
! geom = 'cylrzp'
! geomsetup = 'cylrzp'
!
!--set position of disc patch for coriolis & centrifugal forces
!
 iexternal_force = 5 ! source terms for coriolis and gravity force
!
!--set boundaries
!
 ibound(:) = 3  ! periodic in y,z
 ibound(1) = 5  ! shearing box in x (==R)  
 nbpts = 0      ! use ghosts not fixed
 asize = 1.0    ! box size (corresponds to a in Hawley/Balbus 1992)
 xmin(:) = -0.5*asize  ! set position of boundaries
 xmax(:) = 0.5*asize
 denszero = 1.0
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
 call set_uniform_cartesian(1,psep,xmin,xmax)        ! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*(xmax(1)-xmin(1))*(xmax(1)-xmin(1))
!--CAREFUL WITH (2) totmass = denszero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/float(ntotal) ! average particle mass
!
!--sound speed
!
 przero = 1.e-5
 cs0 = sqrt(gamma*przero/denszero)
 uuzero = przero/(denszero*(gamma-1.))
 polyk0 = przero/denszero**gamma
 write(iprint,*) ' cs0 = ',cs0,' cs2 = ',cs0**2,' u = ',uuzero
 write(iprint,*) ' polyk = ',polyk,' setting to = ',polyk0
 polyk = polyk0

 write(iprint,*) ' Omega_c = ',Omega0, ' R_c = ',Rcentre
 write(iprint,*) ' orbital time = ',2.*pi/Omega0
 write(iprint,*) ' c / (R omega) = ',cs0/(Rcentre*Omega0)

!
!--field strength (beta = pr/(0.5*B^2))
!
 betamhd = 4000.
 Bzeroz = sqrt(2.*przero/betamhd)
 write(iprint,*) ' beta(MHD) = ',betamhd,' Bz = ',Bzeroz,' P = ',przero

!--
 valfvenz = sqrt(Bzeroz**2/denszero)
 wavekmin = sqrt(15.)/4.*Omega0/valfvenz
 
! wavekmin = 2.*pi/((xmax(2)-xmin(2))*Omega0)*sqrt(2.*przero/(denszero*betamhd))
! write(iprint,*) ' qmin = ',wavekmin
 write(iprint,*) 'most unstable wavelength = ',2.*pi/wavekmin

!--setup quantities for the linear inertial waves
 if (imhd.eq.0 .and. isheartest.eq.1) then
    wavenum     = 2.*pi/(xmax(1) -xmin(1))
    wavenumsq   = wavenum*wavenum
    Omega2_wave = 2.*(2.-domegadr)*Omega0*Omega0
    U           = sqrt(cs0*cs0 + Omega2_wave/wavenumsq)
 endif

!
!--now assign particle properties
! 
 iseed = -213
 do i=1,ntotal
    if (imhd.eq.0 .and. isheartest.eq.1) then
       ksi     = ksi0*sin(wavenum*x(1,i))  !--Here, x(1,i) =X0
       dksidX0 = wavenum*ksi0*cos(wavenum*x(1,i))
       x(1,i)  = x(1,i) + ksi !--now x(1,i) = X
    endif
    vel(:,i) = 0. !0.01*(ran1(iseed)-0.5)*cs0
    vel(3,i) = -domegadr*Omega0*x(1,i)
    dens(i) = denszero
    if (imhd.eq.0 .and. isheartest.eq.1) then
       vel(1,i) = vel(1,i) -U*dksidX0
       vel(3,i) = vel(3,i) - (2.-domegadr)*Omega0*ksi
       dens(i)  = dens(i)/(1.+dksidX0)
    endif
    pmass(i) = massp
    uu(i) = uuzero   !*(1. + 0.01*(ran1(iseed)-0.5)) ! isothermal
    if (isheartest.eq.0) then
       select case(ipert)
       case(2)
          vel(1,i) = 0.
       case(1)
          vel(1,i) = 0.01*cs0
       case default
          stop 'invalid ipert'
       end select
    endif
    
    if (imhd.gt.0) then
       Bfield(:,i) = 0.
       Bfield(2,i) = Bzeroz !*sin(2.*pi*x(1,i))
       !if (abs(x(1,i)).lt.0.25*asize) then
       !   Bfield(2,i) = Bzeroz !*(1. + 0.01*sin(4.*2.*pi*(x(2,i)-xmin(2))))
       !else
       !   Bfield(2,i) = 0.
       !endif
       Bfield(3,i) = 0.
    elseif (imhd.lt.0) then
       Bevol(:,i) = 0.
       Bevol(3,i) = 0. !0.5/pi*Bzeroz*cos(2.*pi*x(1,i))
    endif
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
