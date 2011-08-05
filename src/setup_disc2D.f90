!-------------------------------------------------------------------------
! Set up a Keplerian disc with a radial density profile in 2D
!-------------------------------------------------------------------------

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
!
!--define local variables
!            
 implicit none
 integer :: i,npartphi,npartprev
 real :: massp,volume,totmass,pindex,rindex,xi,yi
 real :: denszero,rr,psepold,omegadot,rmin,rmax,rbuffer

 write(iprint,*) 'Accretion disc setup in 2D'
!
!--set co-ordinate system
!
 igeom = 2
 iexternal_force = 2
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 ibound(1) = 2	! reflective in r
 rmin = 0.5	! rmin
 rmax = 1.5  ! rmax
 rbuffer = 0.  ! size of buffer region (gaussian falloff)
 
!
!--set density profile - works for rindex <= -1 and rindex .ne. -2
!
 rindex = -1.5   ! density \propto r**(pindex)
 
 pindex = rindex+2.
 if (ndim.ge.2)then
    ibound(2) = 3  ! periodic in phi
    xmin(2) = 0.0  ! phi_min
    xmax(2) = 2.*pi  ! phi_max
 endif
 if (ndim.ge.3) ibound(3) = 0  ! nothing in z
!
!--set bounds of \eta coordinate according to the index (-1 gives r)
!
 xmin(1) = min((rmin+rbuffer)**pindex,(rmax-rbuffer)**pindex)
 xmax(1) = max((rmin+rbuffer)**pindex,(rmax-rbuffer)**pindex)
 write(iprint,*) 'rmin, rmax = ',rmin,rmax,' coordmin,max = ',xmin(1),xmax(1)
!
!--make sure the particle separation is an even division of the \phi boundary
!
 psepold = psep
 npartphi = INT((xmax(2)-xmin(2))/psep)
 psep = (xmax(2)-xmin(2))/REAL(npartphi)
 write(iprint,*) ' psep = ',psepold,' adjusted to ',psep,' npartphi = ',npartphi
!
!--set up the uniform density grid (uniform in \eta,phi, 
!  where \eta = r**(-(pindex+1))
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--if not cylindrical coordinates, translate to cartesians
!
 do i=1,npart
!
!--use coords in cylindricals to set velocities
!
    rr = x(1,i)**(1./pindex)
    omegadot = SQRT(1./rr**3)
    vel(1,i) = -rr*SIN(x(2,i))*omegadot
    vel(2,i) = rr*COS(x(2,i))*omegadot
    xi = rr*COS(x(2,i))
    yi = rr*SIN(x(2,i))
    x(1,i) = xi
    x(2,i) = yi
!!    CALL coord_transform(x(:,i),ndim,igeom,x(:,i),ndim,1)
 enddo
!
!--then set buffer regions
! 
 npartprev = npart
 print*,'npartprev = ',npartprev
! xmin(1) = 0.
! xmax(1) = LOG((rmin+rbuffer)/rmin)
! call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
! do i=npartprev,npart
!
!--use coords in cylindricals to set velocities
!
!    rr = (rmin+rbuffer)*exp(-x(1,i))
!    omegadot = SQRT(1./rr**3)
!    vel(1,i) = -rr*SIN(x(2,i))*omegadot
!    vel(2,i) = rr*COS(x(2,i))*omegadot
!    xi = rr*COS(x(2,i))
!    yi = rr*SIN(x(2,i))
!    x(1,i) = xi
!    x(2,i) = yi
!!    CALL coord_transform(x(:,i),ndim,igeom,x(:,i),ndim,1)
! enddo

 igeom = 0
 ibound = 0    ! no boundaries

!
!--determine particle mass
!
 denszero = 1.0
 volume = pi*(rmax**2 - rmin**2)
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
!  (note *do not* setup smoothing length as it depends on the conservative
!   density dens which has not yet been calculated) 
!
 do i=1,ntotal
!
!--vel is keplerian
!
    rr = SQRT(DOT_PRODUCT(x(:,i),x(:,i)))
    dens(i) = denszero*(rr**rindex)
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  Exiting subroutine setup'
  
 return
end
