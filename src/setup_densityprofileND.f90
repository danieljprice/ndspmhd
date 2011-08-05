!----------------------------------------------------------------
!     Set up an arbitrary (radial) density profile
!     using random placement of particles
!----------------------------------------------------------------

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
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,idist
 real :: massp,totmass,ran1
 real :: rr,rmass,rsoft
 real, dimension(ndim) :: xpos
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
 npart = int((1./psep)**3)
 print*,'enter npart'
 read*,npart
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
 iseed = -26588
 idist = 2 ! choice of distribution
 select case(idist)
 case(1)
    write(iprint,*) ' Plummer profile'
 case(2)
    write(iprint,*) ' Hernquist model'
 end select
 rsoft = 1.0
 write(iprint,*) ' Total mass = ',totmass,' rsoft = ',rsoft
 
 do i=1,ntotal
!
!--choose random number for mass fraction (0->0.99)
!
    rmass = 0.99*ran1(iseed)
!
!--convert mass fraction to radial position for chosen distribution
!
    call getr_from_m(idist,rmass,rsoft,rr)
    xpos(1) = rr
    if (rr.lt.0.) print*,i,'m = ',rmass,' r = ',rr
!
!--choose random number for (spherical) angle phi (-pi -> pi)
!
    if (ndim.ge.2) xpos(2) = 2.*pi*(ran1(iseed) - 0.5)
!
!--choose random number for (spherical) angle theta (0 -> pi)
!  use ACOS to get even spread, not bunched near poles
!
    if (ndim.ge.3) xpos(3) = ACOS(2.*ran1(iseed) - 1.)
!
!--transform to cartesian coords
!
    call coord_transform(xpos(:),ndim,3,x(:,i),ndim,1)
!
!--set other quantities
!
    vel(:,i) = 0.
    dens(i) = densprofile(idist,rr,rsoft)
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
 
contains

!
! subroutine containing mass distributions to choose from
!
subroutine getr_from_m(ioption,rmass,rsoft,r) 
  implicit none
  integer, intent(in) :: ioption
  real, intent(in) :: rmass,rsoft
  real, intent(out) :: r
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

end subroutine setup
