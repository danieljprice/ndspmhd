!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

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
 integer :: i,iseed
 real :: Rcentre,HoverR,dfac,rhomax,Rtorus
 real :: rmintorus,rmaxtorus,ran1,thetamin,thetamax
 real :: rr,phi,theta,rhoi,rhonorm,przero,fran,zi,rcyl,rfromtorus
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(torus)'
!
!--set boundaries
! 	    
 ibound = 0	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -2.	! set position of boundaries
 xmax(:) = 0.
!
!--set up the uniform density grid
!
 Rcentre = 1.0
 HoverR = 0.1
 dfac = 1.1
 rhomax = 1.0
 Rtorus = Rcentre*HoverR
 
 rmintorus = Rcentre - Rtorus
 rmaxtorus = Rcentre + Rtorus
 thetamin = 0.5*pi - abs(ATAN(0.1))
 thetamax = 0.5*pi + abs(ATAN(0.1))
 
 write(iprint,*) 'Setup for Keplerian torus thing (Owen, Laure)'
 write(iprint,*) ' rmin, max = ',rmintorus,rmaxtorus
 write(iprint,*) ' theta min,max = ',thetamin,thetamax
 write(iprint,*) 'enter particle number '
 read*,npart
 ntotal = npart

 call alloc(ntotal)

 iseed = -7834

 i = 0 
 do while (i < npart)
!
!--choose random number in r, phi and theta
!
   rr = rmintorus + ran1(iseed)*(rmaxtorus-rmintorus)
   phi = 2.*pi*(ran1(iseed) - 0.5)
   theta = thetamin + (ACOS(2.*ran1(iseed) - 1.))/pi*(thetamax-thetamin)
   
   rcyl = rr*COS(theta - 0.5*pi)
   zi = rr*SIN(theta - 0.5*pi)
   rfromtorus = sqrt((rcyl - Rcentre)**2 + zi**2)
   
!
!--compute density from formula
!
   if (rfromtorus.le.Rtorus) then ! check we are in the torus
      rhoi = rhofunc(rr,theta,gamma,rhomax,dfac,Rcentre)
      rhonorm = rhoi/rhomax
   else
      rhonorm = 0.
      rhoi = 0.
   endif
   
   fran = ran1(iseed)
   if (fran < rhonorm) then
      i = i + 1
      x(1,i) = rr*SIN(theta)*COS(phi)
      x(2,i) = rr*SIN(theta)*SIN(phi)
      x(3,i) = rr*COS(theta)
      dens(i) = rhoi
      pmass(i) = 1.0
!      przero = AA*rhoi**gamma
      uu(i) = 0. !! przero/((gamma-1.)*rhoi)
      Bfield(:,i) = 0.
   endif
   
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return

contains

real function rhofunc(rr,theta,gamma,rhomax,dd,R0)
 implicit none
 real, intent(in) :: rr,theta,gamma,rhomax,dd,R0
 real :: polyn,AA,term
!
!--work out A from rhomax
!
 polyn = 1./(gamma-1.)
 AA = 1./((polyn + 1.)*rhomax**(1./polyn))*((dd-1.)/(2.*dd))
!
!--now functional form of rho
!
 term = (1./((polyn + 1.)*R0*AA) &
          *(R0/rr - 0.5*(R0/(rr*SIN(theta)))**2 - 1./(2.*dd)))
 rhofunc = term**polyn

end function rhofunc

end

