!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for a uniform spherical distribution in 1, 2 or 3 dimensions    !!
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
 REAL :: denszero,uuzero,massp

 WRITE(iprint,*) 'Uniform spherical distribution'
!
!--set bounds of initial setup
!                   
 rmax = 1.0
 ibound = 0	! no boundaries 
 iexternal_force = 1	! use toy star force
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
 denszero = totmass/totvol
 massp = totmass/REAL(npart)
 uuzero = 0.1
 WRITE(iprint,10) denszero,uuzero
10 FORMAT(/,' initial density = ',f7.3, ' initial u = ',f7.3,/)
!
!--set these for all particles
! 
 vel(:,:) = 0.
 dens(:) = denszero
 uu(:) = uuzero
 pmass(:) = massp
 Bfield(:,:) = 0.
 WHERE (dens > 0.)
    hh = hfact*(massp/dens(:))**hpower
 END WHERE
 
 RETURN
END SUBROUTINE setup
