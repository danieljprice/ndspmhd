!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for a polytrope in 1, 2 or 3 dimensions                         !!
!!                                                                        !!
!!  Sets up uniform spherical cloud with low temperature                  !!
!!                                                                        !!                                                                    !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
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
 integer :: i,j
 real :: rmax,totmass,totvol
 real :: denszero,uuzero,massp

 write(iprint,*) 'polytrope'
 if (igravity.eq.0) write(iprint,*) 'warning: gravity not set'
!
!--set bounds of initial setup
!                   
 rmax = 2.0
 ibound = 0	! no boundaries 
!
!--setup a uniform sphere of particles
! 
 call set_uniform_spherical(1,rmax)	! 4 = random
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
  
 totmass = totvol
 denszero = totmass/totvol
 massp = totmass/real(npart)
 uuzero = 0.1
 write(iprint,10) denszero,uuzero
10 format(/,' initial density = ',f7.3, ' initial u = ',f7.3,/)
!
!--set these for all particles
! 
 vel(:,:) = 0.
 dens(:) = denszero
 uu(:) = uuzero
 pmass(:) = massp
 bfield(:,:) = 0.
 
 return
end subroutine setup
