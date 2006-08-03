!-------------------------------------------------------------------------
! Set up a particles on rings (useful for discs) with a
! given radial density profile in 2D
!-------------------------------------------------------------------------

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
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer, parameter :: maxrings = 100
 integer :: i,iring,ipart,npartphi,nringsguess,nrings
 real, dimension(maxrings) :: rhoring,drhoring,rring
 real :: massp,volume,totmass,rindex
 real :: denszero,rhoringi,dr,rr,psepold,omegadot,rmin,rmax,rbuffer
 real :: phi,phimin,phimax
 real :: rhor,drhor
 real :: gamm1

 write(iprint,*) 'Accretion disc setup in 2D'

 iexternal_force = 2
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 ibound = 0	! no boundaries
 rmin = 0.5	! rmin
 rmax = 1.5     ! rmax
 rbuffer = 0.3   ! size of buffer region (quadratic falloff)
 phimin = 0.0
 phimax = 2.*pi ! this must be 2\pi unless in cylindrical co-ords
!
!--set density profile
!
 denszero = 1.0
 rindex = -1.5   ! density \propto r**(rindex)
 gamm1 = gamma -1.
 if (gamm1.lt.1.e-3) stop 'error: isothermal eos not implemented'
 write(iprint,*) ' sound speed = ',sqrt(gamma*polyk*denszero**gamm1)
!
!--particle separation determines number of particles in each ring
!
 psepold = psep
 npartphi = INT((phimax-phimin)/psep)
 psep = (phimax-phimin)/REAL(npartphi) ! adjust so exact division of 2\pi
!
!--setup radius of all the rings
!
 iring = 0
 rr = rmin + rbuffer
 do while (rr.lt.rmax)
    iring = iring + 1 
    rring(iring) = rr
    rhoring(iring) = rhor(rindex,rr,rmin,rmax,rbuffer) 
    drhoring(iring) = drhor(rindex,rr,rmin,rmax,rbuffer)
    dr = psep/(rr*rhoring(iring))  ! dr varies depending on the radial density 
    rr = rr + dr                   ! profile (dr=const for 1/r)
 enddo
 nrings = iring
!
!--now go backwards from rmin+rbuffer to rmin
!
 rr = rring(1) - psep/(rring(1)*rhoring(1))
 do while (rr.gt.rmin)
    iring = iring + 1
    rring(iring) = rr
    rhoring(iring) = rhor(rindex,rr,rmin,rmax,rbuffer) 
    drhoring(iring) = drhor(rindex,rr,rmin,rmax,rbuffer)
    dr = psep/(rr*rhoring(iring))
    rr = rr - dr
 enddo 
 nrings = iring
 write(iprint,*) ' number of rings = ',nrings,' particles per ring = ',npartphi 
!
!--allocate memory
!
 call alloc(nrings*npartphi)
!
!--determine particle mass
!
 ntotal = nrings*npartphi
 volume = pi*(rmax**2 - rmin**2)
 totmass = denszero*volume
 massp = totmass/FLOAT(ntotal) ! average particle mass

!---------------------------------------------
!  now setup the particles on these rings
!---------------------------------------------
 do iring=1,nrings
    rr = rring(iring)
    rhoringi = rhoring(iring)
    write(iprint,"(1x,a,i2,1x,2(a,f6.2))") ' ring ',iring,' r = ',rr,' rho = ',rhoringi
    print*,drhoring(iring)
    !
    !--for each ring, set up same number of particles from 0-> 2\pi
    !
    do i=1,npartphi
       ipart = ipart + 1
       npart = npart + 1
       if (ipart.gt.SIZE(dens)) call alloc(npart+npartphi)
       phi = phimin + (i-1)*psep
       !
       !--translate r, phi back to cartesian co-ords
       !
       x(1,ipart) = rr*COS(phi)
       x(2,ipart) = rr*SIN(phi)
       dens(ipart) = rhoringi*denszero
       !
       !--keplerian velocity profile balancing pressure gradient
       !
       omegadot = SQRT(1./rr**3 + polyk/rr*dens(ipart)**gamm1*drhoring(iring))
       vel(1,ipart) = -rr*SIN(phi)*omegadot
       vel(2,ipart) = rr*COS(phi)*omegadot
       uu(ipart) = polyk/gamm1*dens(ipart)**gamm1
       Bfield(:,ipart) = 0.
       pmass(ipart) = massp
    enddo

 enddo

 if (ntotal.ne.npart) stop 'something wrong in setup...npart.ne.ntotal'
 if (ntotal.ne.SIZE(dens)) call alloc(ntotal)   ! trim memory if too much
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  Exiting subroutine setup'
  
 return
end

!--------------------------------------------------------
!  this function specifies the radial density profile
!--------------------------------------------------------

real function rhor(rindex,rr,rmin,rmax,rbuffer)
  implicit none
  real :: rindex,rr,rmin,rmax,rbuffer
  real :: r1,r2,rhor1,drhor1
  real :: aa1,bb1,cc1,aa2,bb2,cc2
!
!--work out quadratic fits for buffer regions
!
  if (abs(rbuffer).gt.1.e-3) then
     r1 = rmin+rbuffer
     rhor1 = 1.0
     drhor1 = rindex*(rhor1)/r1
     call fit_quadratic(rmin,r1,0.0,1.0,drhor1,aa1,bb1,cc1)
  
     r2 = rmax-rbuffer
     rhor1 = 1.0
     drhor1 = rindex*(rhor1)/r2
     call fit_quadratic(rmax,r2,0.0,1.0,drhor1,aa2,bb2,cc2)
  else
     r1 = rmin
     r2 = rmax
  endif
!
!--now setup radial density profile
!
  !!print*,'h ', rr,r1,r2,rmin,rmax
  if (rr.lt.r1) then  ! quadratic falloff in buffer regions
     rhor = (aa1*rr**2 + bb1*rr + cc1)*(r1**rindex)
  elseif (rr.gt.r2) then
     rhor = (aa2*rr**2 + bb2*rr + cc2)*(r2**rindex)
  elseif (rr.ge.rmin .and. rr.le.rmax) then
     rhor = rr**rindex
  else
     rhor = 0.
  endif
  
end function rhor

!------------------------------------------------------------------------
!  this function specifies the derivative of the radial density profile
!------------------------------------------------------------------------

real function drhor(rindex,rr,rmin,rmax,rbuffer)
  implicit none
  real :: rindex,rr,rmin,rmax,rbuffer
  real :: r1,r2,rhor1,drhor1
  real :: aa1,bb1,cc1,aa2,bb2,cc2
!
!--work out quadratic fits for buffer regions
!
  if (abs(rbuffer).gt.1.e-3) then
     r1 = rmin+rbuffer
     rhor1 = 1.0
     drhor1 = rindex*(rhor1)/r1
     call fit_quadratic(rmin,r1,0.0,1.0,drhor1,aa1,bb1,cc1)
  
     r2 = rmax-rbuffer
     rhor1 = 1.0
     drhor1 = rindex*(rhor1)/r2
     call fit_quadratic(rmax,r2,0.0,1.0,drhor1,aa2,bb2,cc2)
  else
     r1 = rmin
     r2 = rmax
  endif
!
!--now setup radial density profile
!
  !!print*,'h ', rr,r1,r2,rmin,rmax
  if (rr.lt.r1) then  ! quadratic falloff in buffer regions
     drhor = (2.*aa1*rr + bb1)*(r1**rindex)
  elseif (rr.gt.r2) then
     drhor = (2.*aa2*rr + bb2)*(r2**rindex)
  elseif (rr.ge.rmin .and. rr.le.rmax) then
     drhor = rindex*rr**(rindex-1)
  else
     drhor = 0.
  endif
  
end function drhor

!-------------------------------------------------------------------
! small utility to fit a quadratic function 
! given two points (x0,y0), (x1,y1) and a gradient dydx(x1)
! returns the coefficients of y = ax^2 + bx + c
!-------------------------------------------------------------------
subroutine fit_quadratic(x0,x1,y0,y1,dydx1,aa,bb,cc)
 implicit none
 real, intent(in) :: x0,x1,y0,y1,dydx1
 real, intent(out) :: aa,bb,cc

 aa = ((y1 - y0) - dydx1*(x1-x0))/(x1**2 - x0**2 - 2.*x1*(x1-x0))
 bb = dydx1 - 2.*aa*x1
 cc = y0 - aa*x0**2 - bb*x0

!! print 10,aa,bb,cc
10 format(/,'quadratic: y = ',f9.3,' x^2 + ',f9.3,' x + ',f9.3)

 return
end subroutine fit_quadratic
