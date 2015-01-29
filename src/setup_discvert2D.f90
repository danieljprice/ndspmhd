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

!
!--determine particle mass
!
 denszero = 1.e-6   ! midplane density
 rhomin   = 1.e-12  ! minimum density
 partvol = psep**ndim
 massp = partvol*denszero ! average particle mass
 omega = sqrt(Mstar/Rdisc**3)
 !H0 = 0.3921651256
 print*,' Rdisc = ',Rdisc,' Mstar = ',Mstar
 
 cs0 = H0*omega !0.0453 !H0*omega
 polyk = cs0**2
 print*,' cs0 = ',cs0
 H0 = cs0/omega
 !cs_physical = cs0*
 !print*,' temperature = ',1.38e-13*,' assuming distance in AU, mass in Msun and G=1'
 print*,' polyk = ',polyk
 print*,' scale height H0 = ',H0

 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.,stretchdim=2,stretchfunc=rhoy) 
 else
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 endif
 npart = ntotal
 print*,'npart =',npart
!
!--set particle mass from actual integral of surface density profile
!
 totmass = 2.*denszero*sqrt(0.5*pi)*H0*erf(xmax(2)/(sqrt(2.)*H0))
 print*,' totmass = ',totmass, ' massp = ',totmass/npart,massp
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
