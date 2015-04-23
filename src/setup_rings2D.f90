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
 real :: massp,volume,totmass,massring
 real :: denszero,rhoringi,dr,rr,psepold,rmin,rmax
 real :: phi,phimin,phimax
 real :: rhor,drhor
 real :: gamm1

 write(iprint,*) 'rings setup in 2D'
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 ibound = 0	! no boundaries
 rmin = 0.0	! rmin
 rmax = 1.0     ! rmax
 phimin = 0.0
 phimax = 2.*pi ! this must be 2\pi unless in cylindrical co-ords
!
!--set density profile
!
 denszero = 1.0
 gamm1 = gamma -1.
 if (gamm1.lt.1.e-3) stop 'error: isothermal eos not implemented'
 write(iprint,*) ' sound speed = ',sqrt(gamma*polyk*denszero**gamm1)
!
!--particle separation determines number of particles in each ring
!
 psepold = psep
 npartphi = INT((phimax-phimin)/psep)
 psep = (phimax-phimin)/REAL(npartphi) ! adjust so exact division of 2\pi

 nrings = 10

 volume = pi*(rmax**2 - rmin**2)
 totmass = denszero*volume
 massring = totmass/nrings
!
!--setup radius of all the rings
!  start at r=midway between rmin and rmax
!
 iring = 0
 rr = 0.5*(rmin + rmax)
!
!--work out rings towards centre
!
 do while (rr.gt.rmin)
    iring = iring + 1 
    rring(iring) = rr
    rhoring(iring) = rhor(rr,rmin,rmax) 
    drhoring(iring) = drhor(rr,rmin,rmax)
    dr = psep/(rr*rhoring(iring))  ! dr varies depending on the radial density 
    rr = rr - dr                   ! profile (dr=const for 1/r)
    print*,iring,'r = ',rr,' dr = ',-dr
 enddo
 nrings = iring - 1
 iring = iring - 1
 print*,'inward rings = ',nrings
!
!--work out rings going outwards
!
 rr = 0.5*(rmin + rmax)
 dr = psep/(rr*rhoring(1))
 rr = rr + dr
 
 do while (rr.lt.rmax)
    iring = iring + 1 
    rring(iring) = rr
    rhoring(iring) = rhor(rr,rmin,rmax) 
    drhoring(iring) = drhor(rr,rmin,rmax)
    dr = psep/(rr*rhoring(iring))  ! dr varies depending on the radial density 
    rr = rr + dr                   ! profile (dr=const for 1/r)
    print*,iring,'r = ',rr,' dr = ',dr
 enddo
 nrings = iring - 1

 write(iprint,*) ' number of rings = ',nrings,' particles per ring = ',npartphi 
!
!--allocate memory
!
 call alloc(nrings*npartphi)
!
!--determine particle mass
!
 ntotal = nrings*npartphi
 massp = totmass/FLOAT(ntotal) ! average particle mass

!---------------------------------------------
!  now setup the particles on these rings
!---------------------------------------------
 ipart = 0
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
       pmass(ipart) = massp
       vel(:,ipart) = 0.
       Bfield(:,ipart) = 0.
    enddo

 enddo

 if (ntotal.ne.npart) stop 'something wrong in setup...npart.ne.ntotal'
 print*,'ntotal = ',ntotal
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

real function rhor(rr,rmin,rmax)
  implicit none
  real :: rr,rmin,rmax

  if (rr.ge.rmin .and. rr.le.rmax) then
     rhor = 1. - rr**2
  else
     rhor = 0.
  endif
  
end function rhor

!------------------------------------------------------------------------
!  this function specifies the derivative of the radial density profile
!------------------------------------------------------------------------

real function drhor(rr,rmin,rmax)
  implicit none
  real :: rr,rmin,rmax
  
  if (rr.ge.rmin .and. rr.le.rmax) then
     drhor = - 2.*rr
  else
     drhor = 0.
  endif  
  
end function drhor

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
