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

!------------------------------------------
! Set up Choi damping test
! Choi, Kim & Witta (2009), ApJS 181, 413
!------------------------------------------

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
 use eos
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,ny
 real :: massp,volume,totmass
 real :: denszero,ampl,vA
 real, dimension(3) :: Bzero
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(choi)'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(:) = 1.
 ny = 6
 if (ndim >= 2) xmax(2) = xmin(2) + ny*psep
 if (ndim >= 3) xmax(3) = xmin(3) + ny*psep
 gamma = 1.
 ampl = 0.1
 gamma_ambipolar = 1000.
 rho_ion = 0.1
!
!--set up the uniform density grid
!
 Bzero(1) = 1.
 Bzero(2) = 0.
 Bzero(3) = 0.
 denszero = 1.0
 vA = sqrt(dot_product(Bzero,Bzero)/denszero)

 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do i=1,npart
    vel(:,i) = 0.
    vel(3,i) = ampl*vA*sin(2.*pi*(x(1,i)-xmin(1)))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0
    Bfield(:,i) = Bzero
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' exiting subroutine setup'
  
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
