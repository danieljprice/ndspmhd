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

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE eos
 USE options
 USE part
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j
 REAL :: rmax,totmass,totvol
 REAL :: rhozero,uuzero,massp

 WRITE(iprint,*) 'Polytrope'
 IF (igravity.EQ.0) WRITE(iprint,*) 'WARNING: gravity not set'
!
!--set bounds of initial setup
!                   
 rmax = 2.0
 ibound = 0	! no boundaries 
!
!--setup a uniform sphere of particles
! 
 CALL set_uniform_spherical(1,rmax)	! 4 = random
!
!--set particle properties
! 
 SELECT CASE (ndim)
  CASE(1)
    totvol = 2.*rmax
  CASE(2)
    totvol = pi*rmax**2
  CASE(3)
    totvol = 4./3.*pi*rmax**3
 END SELECT
  
 totmass = totvol
 rhozero = totmass/totvol
 massp = totmass/REAL(npart)
 uuzero = 0.1
 WRITE(iprint,10) rhozero,uuzero
10 FORMAT(/,' initial density = ',f7.3, ' initial u = ',f7.3,/)
!
!--set these for all particles
! 
 vel(:,:) = 0.
 rho(:) = rhozero
 uu(:) = uuzero
 enin(:) = uuzero
 pmass(:) = massp
 Bfield(:,:) = 0.
 WHERE (rho > 0.)
    hh = hfact*(massp/rho(:))**hpower
 END WHERE
 
 RETURN
END SUBROUTINE setup
