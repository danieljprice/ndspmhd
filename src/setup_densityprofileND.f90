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
!     Set up an arbitrary (radial) density profile
!     using random placement of particles
!----------------------------------------------------------------
module plummer_setup
 implicit none
 real :: mass1,mass2,rsoft,rsoft2,xdist2
 integer :: npart1,npart2
 integer, parameter :: idist = 2 
 
end module plummer_setup


subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 use options, only:ibound
 use part
 use setup_params, only:pi,psep
 use geometry, only:coord_transform
 use timestep, only:iseedMC
 use mem_allocation, only:alloc
 use plummer_setup
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,j
 real :: massp,totmass,ran2,massratio
 real :: rr,rmass
 real, dimension(ndim) :: xpos
 real :: dxij,dxmean,q,vmax,vesc,vr,rr2
 logical, parameter :: isetvelocities = .false.
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(density profile)'
!
!--set boundaries
! 	    
 ibound = 0
 nbpts = 0

 write(iprint,*) ' Radial density profile'
!
!--set up the uniform density grid
!
 npart = 100000
! npart = int((1./psep)**3)
! print*,'enter npart'
! read*,npart
 call alloc(npart)
 
 ntotal = npart
 write(iprint,*) 'npart =',npart
!
!--determine particle mass
!
 totmass = 1.0
 massp = totmass/float(npart) ! average particle mass
!
!--now assign particle properties
! 
!--initialise random number generator
! iseed = -26588
 iseed = iseedMC
 write(iprint,*) ' iseed = ',iseed
 !!idist = 2 ! choice of distribution
 select case(idist)
 case(1)
    write(iprint,*) ' Plummer profile'
 case(2)
    write(iprint,*) ' Hernquist model'
 end select
 rsoft = 1.0
 write(iprint,*) ' Total mass = ',totmass,' rsoft = ',rsoft
 call getr_from_m(idist,0.99,rsoft,rr)
 write(iprint,*) ' Cutoff radius = ',rr

!
!--option for a satellite halo with different parameters
!
 massratio = 1.0  ! 1.0 for no second component, 0.0 for satellite only
 npart1 = int(massratio*npart)
 npart2 = npart - npart1
 mass1 = massratio*totmass
 mass2 = totmass*(1.-massratio)
 if (npart1.ne.npart) then
    rsoft2 = 0.1   ! scale length for 2nd component
    xdist2 = 10.0  ! distance to 2nd component (in x)
    write(iprint,*) ' + smaller satellite, mass ratio = ',massratio
    write(iprint,*) ' mass (1) = ',mass1,' mass(2) = ',mass2
    write(iprint,*) ' npart (1) = ',npart1,' npart(2) = ',npart2
    write(iprint,*) ' rsoft (2) = ',rsoft2,' separation = ',xdist2
    call getr_from_m(idist,0.99,rsoft2,rr)
    write(iprint,*) ' Cutoff radius (2) = ',rr
 endif

 
 do i=1,ntotal
!
!--choose random number for mass fraction (0->0.99)
!
    rmass = 0.99*ran2(iseed)
!
!--convert mass fraction to radial position for chosen distribution
!
    if (i.le.npart1) then
       call getr_from_m(idist,rmass,rsoft,rr)
    else
       call getr_from_m(idist,rmass,rsoft2,rr)
    endif
    xpos(1) = rr
    if (rr.lt.0.) print*,i,'m = ',rmass,' r = ',rr
!
!--choose random number for (spherical) angle phi (-pi -> pi)
!
    if (ndim.ge.2) xpos(2) = 2.*pi*(ran2(iseed) - 0.5)
!
!--choose random number for (spherical) angle theta (0 -> pi)
!  use ACOS to get even spread, not bunched near poles
!
    if (ndim.ge.3) xpos(3) = ACOS(2.*ran2(iseed) - 1.)
!
!--transform to cartesian coords
!
    call coord_transform(xpos(:),ndim,3,x(:,i),ndim,1)
!
!--set other quantities
!
    if (i.le.npart1) then
       dens(i) = densprofile(idist,rr,rsoft)
    else
       !--offset 2nd component by specified amount (in x only)
       x(1,i) = x(1,i) + xdist2 
       dens(i) = densprofile(idist,rr,rsoft2)
    endif
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    Bfield(:,i) = 0.
 enddo
 
 if (isetvelocities) then
    write(iprint,*) 'setting particle velocities...'
    if (npart1.ne.npart) stop 'not implemented for multiple components'
    iseed = -7854
    do i=1,ntotal
   !
   !--set velocities following Aarseth, Henon & Wielen (1974)
   !
       q = 0.0
       vmax = 0.1

       do while (vmax.gt.distfn(q,idist))
          q = ran2(iseed)
          vmax = 0.1*ran2(iseed)
          !!print*,i
       enddo
       rr2 = dot_product(x(:,i),x(:,i))
       vesc = sqrt(2.)*(1. + rr2)**(-0.25)
       vr = q*vesc

       vel(1,i) = (1. - 2.*ran2(iseed))*vr
       q = ran2(iseed)
       vel(2,i) = sqrt(vr**2 - vel(1,i)**2)*COS(2.*pi*q)
       vel(3,i) = sqrt(vr**2 - vel(1,i)**2)*SIN(2.*pi*q)
   !
   !--set velocities to some fraction of equilibrium value
   !
   !    vel(:,i) = 1.2*vel(:,i)
    enddo
 endif
!
!--find mean particle separation
! 
! dxmean = 0.
! do i=1,ntotal
!    do j=i+1,ntotal
!       dxij = sqrt(dot_product(x(:,i)-x(:,j),x(:,i)-x(:,j)))
!       dxmean = dxmean + dxij
!    enddo
! enddo
! dxmean = dxmean/real((ntotal**2 - ntotal)/2)
! write(iprint,*) 'mean particle spacing = ',dxmean
! write(iprint,*) 'suggested softening h = ',dxmean/40.,' to ',dxmean/35.

 if (isetvelocities) then
    call reset_centre_of_mass(x,pmass)
 endif
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
 
contains

!
! subroutine containing mass distributions to choose from
!
subroutine getr_from_m(ioption,rmass,rsoft,rr) 
  implicit none
  integer, intent(in) :: ioption
  real, intent(in) :: rmass,rsoft
  real, intent(out) :: rr
  real :: rmass23
!
!--transformation between a point in the mass distribution and a radial position
!  ie. having solved M(r) = 4\pi \int \rho(r) r^2 dr for r.
!
  select case(ioption)
  case(1)
!
!--Plummer sphere
!
   rmass23 = rmass**(2./3.)
   rr = sqrt(rmass23*rsoft**2/(1. - rmass23))
  case(2)
!
!--Hernquist model
!
   rr = rsoft*(rmass + sqrt(rmass))/(1.-rmass)
  case default
   rr = 0.
   print*,'ERROR in setup, wrong distribution'
  end select
  
end subroutine getr_from_m

!
! Function containing the density profile
!
real function densprofile(ioption,r,rsoft)
  implicit none
  integer :: ioption
  real :: r,rsoft
  
  select case(ioption)
  case(1)
!
!--Plummer sphere
!
   densprofile = 3.*rsoft**2/(4.*pi*(r**2 + rsoft**2)**2.5)

  case(2)
!
!--Hernquist model
!
   densprofile = rsoft/(2.*pi*r*(rsoft + r)**3)
   
  case default
   densprofile = 0.
   print*,'ERROR in setup, wrong distribution'
  
  end select

end function densprofile

!
! Function containing the distribution function for velocities
!
real function distfn(q,ioption)
  implicit none
  integer :: ioption
  real :: q
  
  select case(ioption)
  case(1)
!
!--Plummer sphere
!
   distfn = q**2*(1.-q**2)**3.5
   
  case(2)
!
!--Hernquist model
!
   distfn = 0.
   stop 'ERROR: distribution function not implemented for hernquist profile'
   
  case default
   distfn = 0.
   print*,'ERROR in setup, wrong distribution'
  
  end select

end function distfn

end subroutine setup

!
!--gives solution for all 3 force components for two plummer spheres at differing locations
!
subroutine force_error_densityprofiles(iprofile,npart,xpts,force,residual,errL2,errL1,ierr)
  use plummer_setup
  implicit none
  real, parameter :: pi = 3.1415926536
  integer, intent(in) :: iprofile,npart
  real, intent(in), dimension(3,npart) :: xpts,force
  real, intent(out), dimension(npart) :: residual
  real, intent(out) :: errL2,errL1
  integer, intent(out) :: ierr
  integer :: i
  real, dimension(3) :: x1,x2,dx1,dx2,df,fexacti
  real :: err2,err1,fmax
!
! check for errors
!
  ierr = 0
  if (rsoft.lt.0.) then
     print*,'error: rsoft < 0 in exact_densityprofile'
     ierr = 3
     return
  endif
  
  errL1 = 0.
  errL2 = 0.
  fmax = 0.
  print*,'npart = ',npart,'masses = ',mass1,mass2

  select case(iprofile)  
  case(1)
!
!--Plummer sphere
!
    !--position of 1st component
    x1 = (/0.0,0.0,0.0/)
    !--position of 2nd component
    x2 = (/0.0,0.0,0.0/)
    x2(1) = xdist2

!--potential
!       do i=1,npart
!          yplot(1,i) = -mass1/sqrt(rsoft**2 + dot_product(dx1,dx1)**2) &
!                       -mass2/sqrt(rsoft2**2 + dot_product(dx2,dx2)**2)
!       enddo
!--force
    do i=1,npart
       dx1(:) = xpts(:,i) - x1(:)
       dx2(:) = xpts(:,i) - x2(:)
       fexacti(1) = -mass1*dx1(1)/ &
                    ((rsoft**2 + dot_product(dx1,dx1))**1.5) &
                  - mass2*dx2(1)/ &
                    ((rsoft2**2 + dot_product(dx2,dx2))**1.5)
       fexacti(2) =-mass1*dx1(2)/ &
                    ((rsoft**2 + dot_product(dx1,dx1))**1.5) &
                  - mass2*dx2(2)/ &
                    ((rsoft2**2 + dot_product(dx2,dx2))**1.5)
       fexacti(3) = -mass1*dx1(3)/ &
                    ((rsoft**2 + dot_product(dx1,dx1))**1.5) &
                  - mass2*dx2(3)/ &
                    ((rsoft2**2 + dot_product(dx2,dx2))**1.5)
       df = fexacti(:) - force(:,i)
!       if (i.lt.10) then
!          print*,'force = ',force(:,i), ' fexact = ',fexacti(:)
!       endif
       fmax = max(fmax,maxval(abs(fexacti(:))))
       err2 = dot_product(df,df)
       err1 = sqrt(err2)
       residual(i) = err1/sqrt(dot_product(fexacti(:),fexacti(:)))
       errL2 = errL2 + err2
       errL1 = errL1 + err1
    enddo

    !
    !--normalise errors (use maximum f value)
    !
    print*,'fmax = ',fmax
    if (fmax.gt.tiny(fmax)) then
       errL1 = errL1/(npart) !!*fmax)
       errL2 = errL2/(npart) !!*fmax**2))
    else
       stop 'error normalising errors'
    endif

  case default
     ierr = 1  
  end select
    
  return
end subroutine force_error_densityprofiles

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
