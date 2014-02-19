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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for a uniform spherical distribution in 1, 2 or 3 dimensions    !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 real :: rmax,totmass,totvol
 real :: denszero,uuzero,massp

 write(iprint,*) 'uniform spherical distribution'
!
!--set bounds of initial setup
!                   
 rmax = 1.0
 ibound = 0 ! no boundaries 
!
!--setup a uniform sphere of particles
! 
! call set_uniform_spherical(1,rmax,centred=.true.,perturb=0.2) ! 4 = random
 call set_uniform_spherical(2,rmax,centred=.true.) ! 4 = random
!
!--set particle properties
! 
 select case (ndim)
  case(1)
    totvol = 2.*rmax
  case(2)
    totvol = pi*rmax**2
  case(3)
    totvol = 4./3.*pi*rmax**3
 end select
  
 totmass = 1.0
 denszero = totmass/totvol
 massp = totmass/real(npart)
 uuzero = 0.01
 write(iprint,10) denszero,uuzero
10 format(/,' initial density = ',f7.3, ' initial u = ',f7.3,/)
!
!--set these for all particles
! 
! polyk = 0.001
 vel(:,:) = 0.
 dens(:) = denszero
 if (gamma.lt.1.00001) then
    uu(:) = polyk
 else
    uu(:) = polyk*denszero**(gamma-1.0)/(gamma-1.0)
 endif
 pmass(:) = massp
 Bfield(:,:) = 0.
 print*,' free fall time = ',sqrt(3.*pi/(32.*denszero))
 call reset_centre_of_mass(x,pmass)
 
 return
end subroutine setup

subroutine modify_dump
 use dimen_mhd, only:ndim
 use part, only:npart,x,vel
 use timestep, only:time
 implicit none
 integer :: i
 real :: ampl
!
!--reset time
!
 time = 0.
 ampl = 0.1
!
!--apply radial velocity perturbation
!
 do i=1,npart
    vel(1:ndim,i) = ampl*x(:,i)
 enddo
 
 return
end subroutine modify_dump
