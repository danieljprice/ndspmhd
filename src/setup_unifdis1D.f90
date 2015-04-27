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
!     Set up a uniform density grid of particles in 1D
!----------------------------------------------------------------

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE eos
 USE options
 USE part
 USE setup_params
!
!--define local variables
!            
 IMPLICIT NONE
 INTEGER :: i
 INTEGER :: imax
 REAL :: sigma,massp,totmass,pri
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup(unifdis)'
 IF (ndim .GT. 1) STOP 'unifdis not implemented for ndim > 1'
!
!--initially set up a uniform density grid
! 	    
 ibound = 1	! boundary type
 nbpts = 6	! use ghosts not fixed
 xmin(1) = -0.5	! set position of boundaries
 xmax(1) = 0.5
 imax = INT((xmax(1)-xmin(1))/psep)
 
 totmass = 1.0	!4./3.
 massp = totmass/imax	! average particle mass
 sigma = 1./SQRT(2.)
 iexternal_force = 0
!
!--allocate memory here
!
 CALL alloc(imax)
 
 DO i=1,imax
    x(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    vel(:,i) = 0.
    IF (x(1,i).LT.0.) THEN
       dens(i) = 10.0
       pmass(i) = 10.0*massp
    ELSE
       dens(i) = 1.0
       pmass(i) = massp
    ENDIF
    pri = 1.0
    uu(i) = pri/((gamma-1.)*dens(i))
    IF (imhd.NE.0) THEN 
       Bfield(1,i) = 0.
       Bfield(2,i) = sigma*dens(i)
       Bfield(3,i) = 0.
    ENDIF 
!   print*,i,x(i),dens(i),uu(i),pmass(i),x(i)-x(i-1)
 ENDDO
  
 npart = imax
 ntotal = npart
 
 RETURN
END

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
