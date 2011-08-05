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
 use bound, only:xmin,xmax
 use options, only:ibound
 use part
 use setup_params
 use geometry, only:coord_transform
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,idist
 real :: massp,volume,totmass,ran1
 real :: denszero,rr,rmass,rsoft
 real, dimension(ndim) :: xpos
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(density profile)'
!
!--set boundaries
! 	    
 ibound = 0	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5

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
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
!--initialise random number generator
 iseed = -2658
 idist = 1 ! choice of distribution
 rsoft = 1.0
 
 do i=1,ntotal
!
!--choose random number for mass fraction (0->1)
!
    rmass = ran1(iseed)
!
!--convert mass fraction to radial position for chosen distribution
!
    call getr_from_m(idist,rmass,rsoft,rr)
    xpos(1) = rr
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
   rr = 0.
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
   densprofile = 3./(4.*pi*(r**2 + rsoft**2)**2.5)

  case(2)
!
!--Hernquist model
!
   densprofile = 0.
   
  case default
   densprofile = 0.
   print*,'ERROR in setup, wrong distribution'
  
  end select

end function densprofile

end subroutine setup
