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
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,rhomin
 real :: partvol,denszero,HonR,H0,omega,cs0,totmass
 real :: umass,utime,udist,solarm,au,gg,years,sigma
 real :: t_cour,t_stop,graindens,grainsize,mm,micron
 logical :: equalmass
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 
 if (ndim /= 2) stop 'setup implemented for 2D only'
 if (iexternal_force /= 15) stop 'need iexternal_force=15 for this problem'
!
!--set boundaries
!
 ibound = 0
 ibound(1) = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 gamma = 1.
 HonR = 0.05
 H0 = HonR*Rdisc
 xmin(1) = -0.5  ! set position of boundaries
 xmax(1) = 0.5
 xmin(2) = -3.*H0   ! set position of boundaries
 xmax(2) = 3.*H0
 equalmass = .true.

! constants to translate to physical units
 mm = 0.1       ! 1 mm in cm
 micron = 1.e-4 ! 1 micron in cm
 au = 1.496d13
 solarm = 1.9891d33
 gg = 6.672d-8
 years = 3.1556926d7
 umass = solarm
 udist = 10.*au
 utime = sqrt(udist**3/(gg*umass))
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
 totmass = 2.*denszero*sqrt(0.5*pi)*H0*erf(xmax(2)/(sqrt(2.)*H0))*(xmax(1) - xmin(1))
 print*,' totmass       =',totmass, ' in MSun =',totmass*umass/solarm
 
 ! various information about grain sizes and the stopping time
 sigma = 2.*denszero*sqrt(0.5*pi)*H0
 print*,' Sigma         =',sigma,' in g/cm^2 =',sigma*(umass/udist**2)
 t_cour = 1.2*psep/cs0
 print*,' t_courant     =',t_cour,' t_cour/t_orb =',t_cour/(2.*pi/omega)
 graindens = 3./(umass/udist**3)
 grainsize = 1.*mm/udist
 print*,' grain size    =',grainsize,' in cm     =',grainsize*udist,' in m =',grainsize*udist*0.01
 print*,' grain density =',graindens,' in g/cm^3 =',graindens*(umass/udist**3)

 t_stop = graindens*grainsize/(denszero*cs0)*sqrt(pi/8.)
 print*,' ts (z=0)      =',t_stop,' ts*omega  =',t_stop*omega,' ts/t_cour =',t_stop/t_cour
 t_stop = graindens*grainsize/(denszero*exp(-0.5*(2.)**2)*cs0)*sqrt(pi/8.)
 print*,' ts (z=2H)     =',t_stop,' ts*omega  =',t_stop*omega,' ts/t_cour =',t_stop/t_cour
 print*,' ts*omega (est)=',graindens*grainsize/sigma

 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.,stretchdim=2,stretchfunc=rhoy) 
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
       dens(i) = denszero*exp(-0.5*(x(2,i)/H0)**2)
       pmass(i) = massp
    else
       dens(i) = denszero
       pmass(i) = massp*exp(-0.5*(x(2,i)/H0)**2) + rhomin*partvol
       !if (i.eq.801) print*,pmass(i),i,exp(-0.5*(x(2,i)/H0)**2)
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
 implicit none
 integer :: i,j
 real :: dust_to_gas_ratio
 
 time = 0.
 dust_to_gas_ratio = 1.e-2
 
 if (idust==2) then
    write(iprint,"(/,a,/)") ' > MODIFYING dump by adding dust (two-fluid) <'
    call alloc(2*ntotal)
    do i=1,npart
       itype(i) = itypegas
       j = npart+i
       call copy_particle(j,i)
       x(:,j) = x(:,i)
       vel(:,j) = vel(:,i)
       pmass(j) = pmass(i)*dust_to_gas_ratio
       dens(j)  = dens(i)*dust_to_gas_ratio
       uu(j)    = 0.
       itype(j) = itypedust
    enddo
    npart = 2*npart
 elseif (idust==1 .or. idust==3 .or. idust==4) then
    write(iprint,"(/,a,/)") ' > MODIFYING dump by adding dust (one-fluid) <'
    do i=1,npart
       dustfrac(i) = dust_to_gas_ratio/(1. + dust_to_gas_ratio)
       pmass(i)    = pmass(i)*(1. + dust_to_gas_ratio)
       deltav(:,i) = 0.
    enddo
 endif

end subroutine modify_dump
