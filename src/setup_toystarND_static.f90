!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the toy star static solution in 1, 2 or 3 dimensions        !!
!!                                                                        !!
!!  Sets up a uniform spherical distribution in 1, 2 or 3 dimensions,     !!
!!  which should then be damped to the static toy star solution           !!
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
 USE polyconst
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j
 REAL :: rmax,totmass,totvol
 REAL :: rhozero,uuzero,massp,rhocentre

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
  
 rhocentre = 1.0		! toy star central density 
 totmass = 2*pi*polyk*rhocentre**(gamma-1.)  ! mass of the toy star 
 
 rhozero = totmass/totvol	! initial density
 massp = totmass/REAL(npart)
 uuzero = 0.1
 WRITE(iprint,10) rhocentre,totmass
10 FORMAT(/,' Toy star static solution ',/, &
            '     central density: ',f7.3, ' total mass = ',f7.3,/)
!
!--set these for all particles
! 
 vel(:,:) = 0.
 rho(:) = rhozero
 uu(:) = uuzero
 pmass(:) = massp
 Bfield(:,:) = 0.
 WHERE (rho > 0.)
    hh = hfact*(massp/rho(:))**hpower
 END WHERE
 
 RETURN
END SUBROUTINE setup
