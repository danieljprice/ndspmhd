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
!     Set up an Alfven wave in 2 or 3 dimensions
!     density is constant so this is easy
!     >> should be compiled with ndimB = 3
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
! real, parameter :: pi = 3.1415926536
 real, dimension(3) :: rvec
 real, dimension(ndim) :: runit
 real :: massp,totmass,denszero,gam1,uuzero,przero
 real :: ampl,wk,xlambda,rmax
 real :: valfven
 real :: vparallel,vperp,vz,vperp0,vz0
 real :: bparallel,bperp,bz,bperp0,bz0
 real :: perturb_sin, perturb_cos, Aperp, Az
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup'
!
!--set direction of wave propagation (runit is unit vector in this direction)
!
 rvec(1) = 1.0
 rvec(2) = 0.
 rvec(3) = 0.
 runit(1:ndim) = rvec(1:ndim) 

 write(iprint,*) ' runit = ',runit
!
!--set boundaries
!             
 ibound = 3     ! periodic boundaries
 nbpts = 0      ! no fixed particles
 xmin(:) = 0.0        ! set position of boundaries
 xmax(:) = 1.0
 if (ndim.ge.2) xmax(2:ndim) = 8.*psep                !/runit(:)
 print*,'xmin,xmax = ',xmin,xmax
!
!--read/set wave parameters
! 
 rmax = dot_product((xmax(:)-xmin(:)),runit)
 ampl = 0.001
! write (*,*) 'enter amplitude of disturbance'
! read (*,*) ampl
 
 xlambda = 1.0        !/cos(anglexy)        !*rmax
! write (*,*) 'enter wavelength lambda'
! read (*,*) xlambda
    
 wk = 2.0*pi/xlambda        !         wave number

!
!--setup parameters
!
 vperp0 = 0.1
 vparallel = 0.0
 vz0 = 0.1
 denszero = 1.0
 przero = 0.1
 bparallel = 1.0
 bperp0 = 0.1 
 bz0 = 0.1
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)
 polyk = przero/(denszero**gamma)        ! override setting in input
!
!--initially set up a uniform density grid (also determines npart)
!
 print*,' setting up uniform density grid'
 call set_uniform_cartesian(2,psep,xmin,xmax,.false.)	! 2 = close packed
!
!--determine particle mass
!
 totmass = denszero*product(xmax(:)-xmin(:))
 massp = totmass/float(npart) ! average particle mass
 print*,'npart,massp = ',npart,massp
 
 do i=1,npart        
    perturb_sin = sin(wk*dot_product(x(:,i),runit))
    perturb_cos = cos(wk*dot_product(x(:,i),runit))
    vperp = vperp0*perturb_sin
    vz = vz0*perturb_cos
    bperp = bperp0*perturb_sin
    bz = bz0*perturb_cos
    
    
    vel(1,i) = vparallel*rvec(1) - vperp*rvec(2)
    vel(2,i) = vparallel*rvec(2) + vperp*rvec(1)
    vel(3,i) = vz
    
    dens(i) = denszero
    pmass(i) = massp
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
    uu(i) = uuzero !+ pri/dens(i)*ampl*sin(wk*ri)        ! if not polytropic

    if (imhd.ne.0) then 
       Bfield(1,i) = bparallel*rvec(1) - bperp*rvec(2)
       Bfield(2,i) = bparallel*rvec(2) + bperp*rvec(1)
       Bfield(3,i) = bz
       if (imhd.lt.0) then
          Aperp = bz0*perturb_sin/(2.*pi)
          Az = bperp0*perturb_cos/(2.*pi)
          Bevol(1,i) = -Aperp*runit(2)
          Bevol(2,i) = Aperp*runit(1)
          Bevol(3,i) = Az
       endif
    endif 
 enddo

 bconst(:) = 0.
 bconst(1) = bparallel*rvec(1)
 bconst(2) = bparallel*rvec(2)

 ntotal = npart
 
 valfven = sqrt(bparallel**2/denszero)

!----------------------------

 write(iprint,*) ' wave set: amplitude = ',ampl,' wavelength = ',xlambda,' k = ',wk
 write(iprint,*) ' alfven speed = ',valfven
 
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
