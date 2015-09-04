!----------------------------------------------------------------
!     Set up an arbitrary density distribution in 1D
!     as a test for higher-dimensional grid->SPH data conversion
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
 use mem_allocation, only:alloc
 use uniform_distributions
 use random, only:ran1
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,ipart,itry
 integer, parameter :: npts = 100
 real :: massp,totmass,rhofunc
 real :: myran,rhomax,xi,dx,rhoi,rhoprev,dxmin
 real, dimension(npts) :: rhogrid
 real, external :: rhointerp,rhoexclude
 logical, external :: iallowed
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(anyrho)'
!
!--set boundaries
! 	    
 ibound = 0	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = 0.	! set position of boundaries
 xmax(:) = 1.

 print*,' enter ntotal '
 read*,ntotal
 call alloc(int(1.1*ntotal))
 
 totmass = 0.
 !--integrate to get total mass
 dx = (xmax(1)-xmin(1))/real(npts)
 rhomax = 0.
 rhoprev = 0.
 
!--setup grid based data
 do i=1,npts
    xi = (i-0.5)*dx
    rhogrid(i) = rhofunc(xi)
    print*,'xgrid = ',xi
 enddo
!--integrate to get total mass
 rhoprev = rhointerp(0.,rhogrid,dx,npts) 
 do i=1,npts
    xi = i*dx
    rhoi = rhointerp(xi,rhogrid,dx,npts)
    totmass = totmass + 0.5*(rhoi + rhoprev)*dx
    rhomax = max(rhomax,rhoi)
    rhoprev = rhoi
 enddo
 print*,'totmass = ',totmass,' rhomax = ',rhomax
!--set particle mass based on total mass
 massp = totmass/ntotal
 
 ipart = 0
 iseed = -233
!--initialise sobol sequence
 call sobseq(-123,xi)
 
 itry = 0
 dxmin = (massp/rhomax)**dndim
 do while (ipart < ntotal)
    itry = itry + 1
    xi = 0.001*(itry-1)*dxmin
    if (xi.gt.xmax(1)) then
       ntotal = ipart
    endif
!    call sobseq(1,xi)
    !xi = ran1(iseed)
    
    rhoi = rhointerp(xi,rhogrid,dx,npts)
    if (iallowed(xi,rhoi,ipart,massp,x(1,1:ipart),hh(1:ipart))) then
       ipart = ipart + 1
       x(1,ipart) = xi
       print*,ipart,xi,rhoi
       hh(ipart) = (massp/rhointerp(xi,rhogrid,dx,npts))
    endif
    !rhoi = rhoexclude(xi,rhogrid,dx,npts,ipart,massp,x(1,1:ipart),hh(1:ipart))
    !myran = ran1(iseed)
    !if (myran < rhoi/rhomax) then
    !   ipart = ipart + 1
    !   x(1,ipart) = xi
    !   print*,ipart,xi,rhoi
    !   hh(ipart) = hfact*(massp/rhointerp(xi,rhogrid,dx,npts))
    !endif
 enddo

 npart = ipart
 ntotal = npart
 print*,'npart =',npart 
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    !!!vel(1,i) = x(1,i)
    dens(i) = rhointerp(x(1,i),rhogrid,dx,npts)
    pmass(i) = massp
    uu(i) = 1.0	! isothermal
    bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end

!
!--get position of particle assuming a linear
!  density profile between two adjacent points
!
subroutine getx(fracm,rho1,rho2,dx,xfrac)
 implicit none
 real, intent(in) :: fracm,rho1,rho2,dx
 real, intent(out) :: xfrac
 real :: AA, BB, CC, drho
 
 drho = (rho2 - rho1)
 if (dx.gt.epsilon(dx)) then
    AA = 0.5*drho/dx
 else
    xfrac = 0.
    return
 endif
 BB = rho1
 if (BB.lt.0.) stop 'rho -ve on input to getx'
 CC = -fracm
 
 if (abs(AA).lt.epsilon(AA)) then ! linear equation
    xfrac = -CC/BB
 else
    xfrac = 0.5/AA*(-BB + sqrt(BB**2 - 4.*AA*CC))
 endif
 return
end subroutine getx

real function rhofunc(xi)
 use setup_params, only:pi
 implicit none
 real, intent(in) :: xi
 
 rhofunc = 2. - xi**2
 
end function rhofunc

real function rhointerp(xi,rhogrid,dxgrid,nx)
 implicit none
 integer, intent(in) :: nx
 real, intent(in), dimension(nx) :: rhogrid
 real, intent(in) :: xi,dxgrid
 integer :: i,ip1
 real :: xgrid,dxfrac

! xi = (i-0.5)*dx
 
 i = int(xi/dxgrid + 0.500001)
 xgrid = (i-0.5)*dxgrid
 dxfrac = (xi - xgrid)/dxgrid
 !print*,'> x = ',xi,xgrid,'cell=',i,xi/dxgrid + 0.5,'frac=',dxfrac,'<'
 
 ip1 = i + 1
 if (ip1.gt.nx) ip1 = ip1 - nx
 if (i.lt.1) i = i+nx
 
 rhointerp = (1.-dxfrac)*rhogrid(i) + dxfrac*rhogrid(ip1)
 
end function rhointerp

!
!--function to get probability EXCLUDING particles
!  we have already set up
!
real function rhoexclude(xi,rhogrid,dxgrid,nx,np,massp,xp,hp)
 implicit none
 integer, intent(in) :: nx
 real, intent(in), dimension(nx) :: rhogrid
 real, intent(in) :: xi,dxgrid
 integer, intent(in) :: np
 real, intent(in) :: massp
 real, intent(in), dimension(np) :: xp,hp
 integer :: i
 real :: rhoin,rhotemp,rhointerp,wabi,h1i,h2i,r2,q2,q,hmin
 real, parameter :: cnormk = 2./3.

 rhoin = rhointerp(xi,rhogrid,dxgrid,nx)
 hmin = 0.5*(massp/rhoin)
 rhotemp = rhoin
 !--subtract off density of particles already set
 overi: do i=1,np
    r2 = (xi - xp(i))**2
    !--exclude if closer than minimum allowed separation
    if (r2.lt.hmin*hmin) then
       rhotemp = 0.
       exit overi
    else
       h1i = 1./hp(i)
       h2i = h1i*h1i
       q2 = r2*h2i
       if (q2 < 4.) then
          q = sqrt(q2)
          if (q2.lt.1.0) then
             wabi = 1. - 1.5*q2 + 0.75*q2*q
          else
             wabi = 0.25*(2.-q)**3
          endif
          rhotemp = max(rhotemp - massp*cnormk*wabi*h1i,0.)
       endif
    endif
 enddo overi
 rhoexclude = rhotemp
 !if (rhoexclude.lt.rhoin) print*,xi,np,' rho = ',rhoin,rhoexclude
 
end function rhoexclude

logical function iallowed(xi,rhoi,np,massp,xp,hp)
 implicit none
 integer, intent(in) :: np
 real, intent(in) :: xi,rhoi,massp
 real, intent(in), dimension(np) :: xp,hp
 integer :: i
 real :: hmin,r2

 hmin = (massp/rhoi)
 iallowed = .true.
 !--exclude if any particle is closer than the allowed minimum separation
 overi: do i=1,np
    r2 = (xi - xp(i))**2
    !--exclude if closer than minimum allowed separation
    if (r2.lt.hmin*hmin) then
       iallowed = .false.
       exit overi
    endif
 enddo overi
 
end function iallowed

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
