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
 real, dimension(ndim) :: runit
 real :: massp,totmass,denszero,gam1,uuzero,przero
 real :: anglexy,ampl,wk,xlambda,rmax
 real :: valfven
 real :: vparallel,vperp,vz,vperp0,vz0
 real :: bparallel,bperp,bz,bperp0,bz0
 real :: perturb_sin, perturb_cos
 real :: Aperp,Az
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup'
!
!--set direction of wave propagation (runit is unit vector in this direction)
!
 anglexy = 30.        ! angle in degrees x,y plane
! anglez = 45.        ! angle in degrees z plane
 anglexy = anglexy*pi/180.        ! convert to radians
 runit(1) = cos(anglexy)
 runit(2) = sin(anglexy)
! runit(3) = 0.
 write(iprint,*) ' runit = ',runit
!
!--set boundaries
!             
 ibound = 3        ! periodic boundaries
 nbpts = 0        ! no fixed particles
 xmin(:) = 0.0        ! set position of boundaries
 xmax(:) = 1.0/runit(:)
 print*,'xmin,xmax = ',xmin,xmax,(xmax(1)-xmin(1))/8
!
!--read/set wave parameters
! 
 rmax = sqrt(dot_product(xmax(:)-xmin(:),xmax(:)-xmin(:)))
 ampl = 0.001
 print*,'rmax = ',rmax
! write (*,*) 'enter amplitude of disturbance'
! read (*,*) ampl
 
 xlambda = 1.0 !!rmax !!!2.0 !!!(xmax(1)-xmin(1))/cos(anglexy)
 write(iprint,*) 'xlambda = ',xlambda
! write (*,*) 'enter wavelength lambda'
! read (*,*) xlambda
    
 wk = 2.0*pi/xlambda        !         wave number


!
!--setup parameters
!
 przero = 0.1

 vperp0 = 0.1
 vparallel = 0.0
 vz0 = 0.1
 denszero = 1.0
 bparallel = 1.0
 bperp0 = 0.1 
 bz0 = 0.1
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)
!
!--initially set up a uniform density grid (also determines npart)
!
 print*,' setting up uniform density grid'
 call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)        ! 2 = close packed
!
!--determine particle mass
!
 totmass = denszero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/float(npart) ! average particle mass
 print*,'npart,massp = ',npart,massp
 
 do i=1,npart        
    perturb_sin = sin(wk*dot_product(x(:,i),runit))
    perturb_cos = cos(wk*dot_product(x(:,i),runit))
    vperp = vperp0*perturb_sin
    vz = vz0*perturb_cos
    bperp = bperp0*perturb_sin
    bz = bz0*perturb_cos
    
    
    vel(1,i) = vparallel*runit(1) - vperp*runit(2)
    vel(2,i) = vparallel*runit(2) + vperp*runit(1)
    vel(3,i) = vz
    
    dens(i) = denszero
    pmass(i) = massp
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
    uu(i) = uuzero !+ pri/dens(i)*ampl*sin(wk*ri)        ! if not polytropic

    if (imhd.ne.0) then
       Bfield(1,i) = Bparallel*runit(1) - Bperp*runit(2)
       Bfield(2,i) = Bparallel*runit(2) + Bperp*runit(1)
       Bfield(3,i) = Bz
       if (imhd.lt.0) then
          Aperp = bz0*perturb_sin/(2.*pi)
          Az = bperp0*perturb_cos/(2.*pi)
          Bevol(1,i) = -Aperp*runit(2)
          Bevol(2,i) = Aperp*runit(1)
          Bevol(3,i) = Az
       endif
    endif 
 enddo
!
!--use constant field subtraction as instability correction
!
 bconst(:) = 0.
 bconst(1) = bparallel*runit(1)
 bconst(2) = bparallel*runit(2)

 ntotal = npart
 
 valfven = sqrt(bparallel**2/denszero)
 if (iener.eq.0) then
    polyk = przero/denszero**gamma
    write(iprint,*) 'setting polyk = ',polyk
 endif


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
