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
!     Set up for an advection test similar to that performed in
!     Stone, Hawley, Evans and Norman (1992) ApJ 388,415-437
!     (see also Evans & Hawley (1988) ApJ 332,659)
!
!     Optional smoothing of discontinuities (not necessary)
!
!     Note for all MHD setups, only the magnetic field should be setup
!----------------------------------------------------------------

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE bound
 USE loguns
 USE eos
 USE options
 USE part
 USE part_in
 USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,ntot
 REAL :: dist,xcentre,vxi,term,term2,denszero,massp,gam1
 REAL :: xpulseleft,xpulseright,deltaleft,deltaright
 REAL :: exx,uuzero,Bypulse,velpulse,denspulse,przero
 REAL :: ampl,spsoundi,velzero,dsmooth
 REAL, DIMENSION(ndimB) :: Bzero

 IF (ndim.NE.1) STOP ' advection test only setup for 1D : recompile '
 
 ibound = 3	! periodic
! nbpts = 10
 xmin(1) = -0.5
 xmax(1) = 0.5
 xcentre  = 0.0
 dsmooth = 0.0		! amount by which to smooth discontinuities ( < 20 )
 xpulseleft = xcentre - 25*psep	! 50 particles wide
 xpulseright = xcentre + 25*psep
 gam1 = gamma - 1.
 dist  = 0.5*psep
 denszero =  10.0   
 massp = denszero*psep                 ! the mass per sph particle
 npart = 0
 uuzero = 1.0/(gamma*gam1)
 ampl = 0.01		! amplitude of the density perturbation
 IF (imhd.NE.0) THEN
    Bypulse = 0.001
 ELSE
    WRITE(iprint,*) 'MHD not set, this test does nothing'
!   CALL quit
    Bypulse = 0.
 ENDIF       
 Bzero(1) = 0.
 Bzero(2) = Bypulse
 Bzero(3) = 0.

 velpulse = 5.0
 velzero = velpulse
 denspulse = denszero*(1.+ampl)
!
!--allocate variables
!
 ntot = NINT((xmax(1)-xmin(1))/psep) + 1
 CALL alloc(ntot)

! i = 1
 x(1,1) = xmin(1) + 0.5*psep
 x(1,2) = xmin(1) + psep + 0.5*psep
 dens(1:2) = denszero
 uu(1:2) = uuzero
 pmass(1:2) = massp
 vel(:,1:2) = 0.
 vel(1,1:2) = velzero
 
 Bfield(:,1:2) = 0.

 i = 1
 DO WHILE (x(1,i).lt.xmax(1))  
    i = i + 1
!   npart = npart + 1
    x(1,i) = xmin(1) + (i-1)*psep + 0.5*psep
    deltaleft = (x(1,i) - xpulseleft)/dist
    deltaright = (x(1,i) - xpulseright)/dist
	
    Bfield(:,i) = 0.		! initially zero, reset later
    vel(:,i) = 0.	 
    vel(1,i) = velzero
    dens(i) = denszero

    IF ((deltaleft.gt.dsmooth).AND.(deltaright.le.-dsmooth)) THEN
       Bfield(2,i) = Bypulse
       vel(1,i) = velpulse
!      dens(i) = denspulse
    ELSEIF (deltaleft.le.-dsmooth) THEN
       Bfield(2,i) = 0.
!      vel(1,i) = velzero
       dens(i) = denszero
    ELSEIF ((deltaleft.gt.-dsmooth).AND.(deltaleft.le.dsmooth)) THEN
       exx = exp(deltaleft)
       Bfield(2,i) = (Bypulse*exx)/(1 + exx)
!      vel(1,i) = (velpulse*exx)/(1+exx)
!      dens(i) = (denszero + denspulse*exx)/(1+exx)
    ELSEIF ((deltaright.gt.-dsmooth).AND.(deltaright.le.dsmooth)) THEN
       exx = exp(deltaright)
       Bfield(2,i) = Bypulse - (Bypulse*exx)/(1 + exx)
!      vel(1,i) = velpulse - (velpulse*exx)/(1+exx)
!      dens(i) = denspulse - ((denspulse-denszero)*exx)/(1+exx)
    ELSE
       Bfield(2,i) = 0.  
       dens(i) = denszero
!      vel(1,i) = velzero
    ENDIF

!   x(1,i) = x(1,i-2) + 2.*massp/dens(i-1)	 
    uu(i) = uuzero
    pmass(i) = massp
 ENDDO
 
 npart = i-1
 ntotal = npart
    
 RETURN
END

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
