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

!----------------------------------------------------------------
!     Set up a disc section in 2D
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use part
 use setup_params
 use externf, only:Rdisc,Mstar
 use eos,     only:gamma,polyk
 
 use uniform_distributions
 use dust,    only:init_drag,grain_size,grain_dens,get_tstop,umass,udist,utime,au,&
                   grain_size_cm,grain_dens_cgs,solarm,years
!
!--define local variables
!            
 implicit none
 integer :: i,ierr,idir!,nz
 real :: massp,rhomin
 real :: partvol,denszero,HonR,H0,omega,cs0,totmass
 real :: sigma!,zmin,zmax,dz,zonh
 real :: t_cour,t_stop,rhogas,rhodust
 logical :: equalmass
 real, parameter :: dust_to_gas_ratio = 1.e-2

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup (discvert)'
 
 select case(ndim)
 case(1)
    if (iexternal_force /= 15) stop 'need iexternal_force=15 for this problem'
 case(2)
    if (iexternal_force /= 15) stop 'need iexternal_force=15 for this problem' 
 case default
    stop 'setup implemented for 1D and 2D only'
 end select
!
!--set boundaries
!
 nbpts = 0      ! use ghosts not fixed
 gamma = 1.
 HonR = 0.05
 H0 = HonR*Rdisc

 if (ndim==2) then
    ibound = 0
    ibound(1) = 3     ! boundaries
    xmin(1) = -0.25  ! set position of boundaries
    xmax(1) = 0.25
    idir = 2
 else ! ndim = 1
    ibound = 0
    idir = 1
 endif
 xmin(idir) = -3.*H0   ! set position of boundaries
 xmax(idir) = 3.*H0
 equalmass = .true.

!
!--determine particle mass
!
 denszero = 1.e-3   ! midplane density
 rhomin   = 1.e-12  ! minimum density
 partvol = psep**ndim
 massp = partvol*denszero ! average particle mass
 omega = sqrt(Mstar/Rdisc**3)
 !H0 = 0.3921651256
 print*,' rho0  = ',denszero,' in g/cm^3 = ',denszero*(umass/udist**3)
 print*,' Rdisc = ',Rdisc,' Mstar = ',Mstar
 print*,' Omega = ',omega,' t_orb = ',2.*pi/omega,' in years = ',2.*pi/omega*utime/years
 print*,' udist = ',udist,' utime = ',utime,' umass = ',umass
 print*,' udist (AU) = ',udist/au,' utime (yrs) = ',utime/years,' umass (solarm) = ',umass/solarm
 print*,' utime/orbits = ',omega/(2.*pi)
 print*,' density units = ',umass/udist**3

 cs0 = H0*omega !0.0453 !H0*omega
 polyk = cs0**2
 print*,' cs0   = ',cs0
 H0 = cs0/omega
 !cs_physical = cs0*
 !print*,' temperature = ',1.38e-13*,' assuming distance in AU, mass in Msun and G=1'
 print*,' polyk = ',polyk
 print*,' scaleheight H =',H0
!
!--set particle mass from actual integral of surface density profile
!
 totmass = 2.*denszero*sqrt(0.5*pi)*H0*erf(xmax(idir)/(sqrt(2.)*H0))
 if (ndim==2) totmass = totmass*(xmax(1) - xmin(1))
 print*,' totmass       =',totmass, ' in MSun =',totmass*umass/solarm
!
!--set grain size from Kdrag in input file
!
 if (idust > 0) then
    grain_size_cm = Kdrag

    ! various information about grain sizes and the stopping time
    sigma = 2.*denszero*sqrt(0.5*pi)*H0
    print*,' Sigma         =',sigma,' in g/cm^2 =',sigma*(umass/udist**2)
    t_cour = 1.2*psep/cs0
    print*,' t_courant     =',t_cour,' t_cour/t_orb =',t_cour/(2.*pi/omega)

    grain_dens = grain_dens_cgs/(umass/udist**3)
    grain_size = grain_size_cm/udist
    print*,' grain size    =',grain_size,' in cm     =',grain_size*udist,' in m =',grain_size*udist*0.01
    print*,' grain density =',grain_dens,' in g/cm^3 =',grain_dens*(umass/udist**3)

    t_stop = grain_dens*grain_size/(denszero*cs0)*sqrt(pi/8.)
    print*,' ts (z=0)      =',t_stop,' ts*omega  =',t_stop*omega,' ts/t_cour =',t_stop/t_cour
    t_stop = grain_dens*grain_size/(denszero*exp(-0.5*(2.)**2)*cs0)*sqrt(pi/8.)
    print*,' ts (z=2H)     =',t_stop,' ts*omega  =',t_stop*omega,' ts/t_cour =',t_stop/t_cour
    print*,' ts*omega (est)=',grain_dens*grain_size/sigma
    print*,' ts            =',grain_dens*grain_size/cs0*sqrt(pi/8.),' / density'

    print "(/,a,f6.3)",' CHECK drag routine: gamma = ',gamma
    call init_drag(ierr,gamma)
    rhogas  = denszero
    rhodust = dust_to_gas_ratio*rhogas
    t_stop  = get_tstop(3,rhogas,0.,cs0,0.) ! assume rhodust = 0. here to match above
    print*,' ts (z=0)            =',t_stop
    print*,' equivalent K (z=0)  = ',rhogas*rhodust/(t_stop*(rhogas + rhodust))

    rhogas  = denszero*exp(-0.5*(2.)**2)
    rhodust = dust_to_gas_ratio*rhogas
    t_stop = get_tstop(3,rhogas,0.,cs0,0.) ! assume rhodust = 0. here to match above
    print*,' ts (z=2H)           =',t_stop
    print*,' equivalent K (z=2H) =',rhogas*rhodust/(t_stop*(rhogas + rhodust))
 endif
 
 ! WRITE cs*ts to file
 !nz = 1000
 !zmax = 5.
 !zmin = -5.
 !dz = (zmax - zmin)/real(nz - 1)
 !open(unit=1,file='rescrit.out',status='replace')
 !print*,' WRITING cs*ts to rescrit.out'
 !write(1,"(a)") '# [    z (AU)   ] [ ts*cs (AU) ] [   z (code units)  ] [ ts*cs (code units) ]'
 !do i=1,nz
 !   zonh = zmin + (i-1)*dz
 !   rhogas = denszero*exp(-0.5*(zonh)**2)
 !   rhodust = dust_to_gas_ratio*rhogas
 !   t_stop = get_tstop(3,rhogas,rhodust,cs0,0.)
 !   write(1,*) zonh*H0*udist/AU,t_stop*cs0*udist/au,zonh*H0,t_stop*cs0
 !enddo
 !close(unit=1)

 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.,stretchdim=idir,stretchfunc=rhoy) 
 else
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 endif
 npart = ntotal
 print*,' npart =',npart
 print*,' massp = ',totmass/npart,massp
 massp = totmass/npart
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    if (equalmass) then
       dens(i) = denszero*exp(-0.5*(x(idir,i)/H0)**2)
       pmass(i) = massp
    else
       dens(i) = denszero
       pmass(i) = massp*exp(-0.5*(x(idir,i)/H0)**2) + rhomin*partvol
    endif
    if (gamma > 1.) then
       uu(i) = cs0**2/(gamma*(gamma-1.)) ! adiabatic
    else
       uu(i) = 1.5*polyk ! isothermal
    endif
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return

contains
 real function rhoy(y)
  real, intent(in) :: y

  rhoy = exp(-0.5*(y/H0)**2)

 end function rhoy

end subroutine setup

subroutine modify_dump
 use loguns,       only:iprint
 use part
 use options,        only:idust
 use timestep,       only:time
 use mem_allocation, only:alloc
 use eos,            only:gamma
 use dust,           only:init_drag
 implicit none
 integer :: i,j,ierr
 real, parameter :: dust_to_gas_ratio = 1.e-2
 
 time = 0.
 if (idust==2) then
    write(iprint,"(/,a,/)") ' > MODIFYING dump by adding dust (two-fluid) <'
    call alloc(2*ntotal)
    j = npart
    do i=1,npart
       itype(i) = itypegas
       j = j + 1
       call copy_particle(j,i)
       x(:,j) = x(:,i)
       vel(:,j) = vel(:,i)
       pmass(j) = pmass(i)*dust_to_gas_ratio
       dens(j)  = dens(i)*dust_to_gas_ratio
       uu(j)    = 0.
       itype(j) = itypedust
    enddo
    npart = j
    ntotal = npart
 elseif (idust==1 .or. idust==3 .or. idust==4) then
    write(iprint,"(/,a,/)") ' > MODIFYING dump by adding dust (one-fluid) <'
    do i=1,npart
       dustfrac(i) = dust_to_gas_ratio/(1. + dust_to_gas_ratio)
       pmass(i)    = pmass(i)*(1. + dust_to_gas_ratio)
       deltav(:,i) = 0.
    enddo
 endif
 
 if (idust /= 0) then
    call init_drag(ierr,gamma)
    if (ierr /= 0) stop 'error initialising drag'
 endif

end subroutine modify_dump
