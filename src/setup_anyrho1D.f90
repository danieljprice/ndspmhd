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
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,ipart,npartadd,iadd
 integer, parameter :: npts = 100
 real :: massp,volume,totmass,rhofunc,xpos,rhopartprev,xoffset
 real :: denszero,myran,rhomax,xi,ran1,dx,rhoi,rhoprev,xprev,totmassprev,frac
 real, dimension(npts) :: rhogrid
 real :: rhointerp
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(anyrho)'
!
!--set boundaries
! 	    
 ibound = 3	! boundaries
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
 
 totmass = 0.
 ipart = 0
 rhoprev = 0.
 xi = 0.
 rhoprev = rhointerp(xi,rhogrid,dx,npts)
 rhopartprev = rhoprev

 xprev = xi
 do i=1,npts
    !xprev = xi
    xi = i*dx
    rhoi = rhointerp(xi,rhogrid,dx,npts)
    totmassprev = totmass
    totmass = totmass + 0.5*(rhoi + rhoprev)*dx
    npartadd = int(totmass/massp)
    if (i.eq.npts) npartadd = nint(totmass/massp)
    !if (i.eq.1) then
    !   npartadd = npartadd - 1 ! already added one particle
    !   totmass = totmass - massp
    !endif
    if (npartadd.gt.0) then
       print*,'npartadd = ',npartadd,' totmass/massp = ',totmass/massp
       rhopartprev = rhointerp(xprev,rhogrid,dx,npts)
       do iadd = 1,npartadd
          frac = iadd*massp
          call getx(frac,rhopartprev,rhoi,(xi-xprev),xoffset)
          print*,'frac = ',frac,rhopartprev,rhoi,xi-xprev,xoffset
          xpos = xprev + xoffset
          ipart = ipart + 1
          x(1,ipart) = xpos
          print*,'xi = ',xpos,xi, 'xprev = ',xprev
       enddo
       totmass = totmass - npartadd*massp
       xprev = xpos
    endif
    rhoprev = rhoi
 enddo

 npart = ipart
 ntotal = npart
 print*,'npart =',npart

 print*,' TEST INTERPOLATION : '
 print*,'rhozero= ',rhointerp(0.,rhogrid,dx,npts),rhointerp(1.,rhogrid,dx,npts)
 print*,' rho (0.5*dx) = ',rhogrid(1),rhointerp(0.5*dx,rhogrid,dx,npts)
 print*,' rho (dx) = ',0.5*(rhogrid(1)+rhogrid(2)),rhointerp(dx,rhogrid,dx,npts)
 
 
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    !!!vel(1,i) = x(1,i)
    dens(i) = rhofunc(x(1,i))
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
 
 rhofunc = 2. + 1.99*sin(10.*pi*xi)
 
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
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
