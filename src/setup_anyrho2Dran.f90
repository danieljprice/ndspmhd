!----------------------------------------------------------------
!     Set up an arbitrary density distribution in 2D
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
 use linklist
!
!--define local variables
!            
 implicit none
 integer :: i,j,iseed,ipart,ix,iy,itry,irholevel,nextra
 integer :: icellx,icelly,icellz,icell,nxpts,nypts,ncellslevel
 integer, parameter :: nx = 100, ny = 100, nrholevels = 5
 integer, dimension(nrholevels) :: ihead
 integer, dimension(nx*ny) :: nextcell,iorder
 real, dimension(nx*ny) :: dxmincell,rhocell
 real :: massp,totmass,rowmass,rowmassprev,rhofunc,fac,drho
 real :: myran,rhomax,xi(ndim),dx,dy,rhoi,rhoprev,dxtry,dxmin
 real :: dxstep,xmincell,ymincell
 real, external :: rhointerp,rhoexclude,ran1,ran2
 logical, external :: iallowed
 logical, parameter :: userandom = .false.,adaptmesh=.true.
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
 rhomax = 0.
 !--integrate to get total mass
 dx = (xmax(1)-xmin(1))/real(nx)
 dy = (xmax(2)-xmin(2))/real(ny)
 rowmassprev = 0.
 do j=1,ny
    xi(2) = j*dy
    rowmass = 0.
    rhoprev = 0.
    do i=1,nx
       xi(1) = i*dx
       rhoi = rhofunc(xi)
       icell = i + (j-1)*nx
       rhocell(icell) = rhofunc(xi)
       rowmass = rowmass + 0.5*(rhoi + rhoprev)*dx
       rhomax = max(rhomax,rhoi)
       rhoprev = rhoi
    enddo
    totmass = totmass + 0.5*(rowmass + rowmassprev)*dy
    rowmassprev = rowmass
 enddo
 print*,'totmass = ',totmass,' rhomax = ',rhomax
!--set particle mass based on total mass
 massp = totmass/ntotal
 
 if (adaptmesh) then
    !--sort grid cells by density
    call indexx(nx*ny,rhocell,iorder)
    !--bin grid cells by density
    !ihead(:) = -1
    !drho = rhomax/nrholevels
    !do j=1,ny
    !   xi(2) = j*dy
    !   do i=1,nx
    !      xi(1) = i*dx
    !      rhoi = rhofunc(xi)
    !      irholevel = int(rhoi/drho) + 1
    !      icell = i + (j-1)*nx
    !      nextcell(icell) = ihead(irholevel)
    !      ihead(irholevel) = icell
    !      dxmincell(icell) = (massp/rhoi)**dndim
    !   enddo
    !enddo
 endif
 
 ipart = 0
 iseed = -233
 if (userandom) then
    fac = 0.85*sqrt(1./0.834)
    call sobseq(-2,xi)
 else
    fac = 1.0*sqrt(1./0.834)
 endif
 ix = 0
 iy = 0
 itry = 0
 dxmin = (massp/rhomax)**dndim
 xi(1) = xmin(1)
 xi(2) = xmin(2)

 if (adaptmesh) then
    do itemp=nx*ny,1,-1
       icell=iorder(itemp)
       !print*,'cell ',icell,' rho = ',rhocell(icell)
       dxstep = 0.005*(massp/rhocell(icell))**dndim

!    do irholevel = nrholevels,1,-1
!       icell = ihead(irholevel)
!       print*,' rho level = ',irholevel,' rho > ',(irholevel-1)*drho
!       ncellslevel = 0
!       do while (icell.ne.-1)
!          ncellslevel = ncellslevel + 1
!          dxstep = 0.005*dxmincell(icell)
          nextra = 2 ! allow overlap between cells
          nxpts = int(dx/dxstep) + 2*nextra
          nypts = int(dy/dxstep) + 2*nextra
          icellz = 1 !(icell - 1)/(nx*ny) + 1
          icelly = (icell - (icellz-1)*nx*ny - 1)/nx + 1
          icellx = icell - (icellz-1)*nx*ny - (icelly-1)*nx
          xmincell = xmin(1) + max((icellx-1)*dx - nextra*dxstep,0.)
          ymincell = xmin(2) + max((icelly-1)*dy - nextra*dxstep,0.)
          !print*,' dx = ',dx,' here = ',dxmincell(icell),dxstep,dx/dxstep
          !print*,' in cell ',icell,' xmin, ymin = ',xmincell,ymincell,' nx,ny = ',nxpts,nypts
          do iy = 1,nypts
             xi(2) = min(ymincell + (iy-1)*dxstep + 0.01*ran1(iseed)*dxstep,xmax(2))
             do ix = 1,nxpts
                xi(1) = min(xmincell + (ix-1)*dxstep + 0.01*ran1(iseed)*dxstep,xmax(1))
                rhoi = rhofunc(xi)
                if (iallowed(xi,rhoi,ipart,massp,x(:,1:ipart),hh(1:ipart),fac)) then
                   itry = 0
                   ipart = ipart + 1
                   x(:,ipart) = xi
                   print*,ipart,xi,rhoi
                   hh(ipart) = (massp/rhoi)**dndim
                endif             
             enddo
          enddo
!          icell = nextcell(icell)
!       enddo
!       print*,' ncells this level = ',ncellslevel,' npart= ',ipart
!       print*,' npart = ',ipart
    enddo
    massp = massp*ntotal/real(ipart)
 else 
 
 do while (ipart < ntotal)
    itry = itry + 1
    if (.not.userandom) then
       ix = ix + 1
       !--adaptive step size in x
       dxtry = dxmin
       !dxtry = 0.01*(massp/rhofunc(xi(1)))**dndim
       xi(1) = xmin(1) + 0.01*(ix-1)*dxtry
       !dx = 0.5*(dxtry + (massp/rhofunc(xi(1)))**dndim)
       !xi(1) = xmin(1) + 0.01*(ix-1)*dx
       !print*,dx,dxmin

       if (xi(1).gt.xmax(1)) then
          iy = iy + 1
          ix = 0
          xi(1) = xmin(1)
          xi(2) = xmin(2) + 0.01*(iy-1)*dxmin
          !print*,xi
          if (xi(2).gt.xmax(2)) ntotal = ipart
       endif
    else
       do idim=1,ndim
          xi(idim) = ran2(iseed)
       enddo
       !call sobseq(ndim,xi)
       if (mod(itry,10000).eq.0) then
          fac = fac*0.9999
          print*,fac
       !else
       !   fac = 0.9*sqrt(1./0.834)
       endif
       !if (itry.gt.10000) print*,itry,fac
    endif    
    rhoi = rhofunc(xi)
    if (iallowed(xi,rhoi,ipart,massp,x(:,1:ipart),hh(1:ipart),fac)) then
       itry = 0
       ipart = ipart + 1
       x(:,ipart) = xi
       print*,ipart,xi,rhoi
       hh(ipart) = fac*(massp/rhoi)**dndim
    endif
    !rhoi = rhofunc(xi) !rhoexclude(xi,rhogrid,dx,npts,ipart,massp,x(1,1:ipart),hh(1:ipart))
    !myran = ran1(iseed)
    !if (myran < rhoi/rhomax) then
    !   ipart = ipart + 1
    !   x(:,ipart) = xi
    !   print*,ipart,xi,rhoi
    !   hh(ipart) = hfact*(massp/rhoi)
    !endif
 enddo
 
 endif

 npart = ipart
 ntotal = npart
 print*,'npart =',npart 
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    !!!vel(1,i) = x(1,i)
    dens(i) = rhofunc(x(:,i))
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
 use dimen_mhd, only:ndim
 use setup_params, only:pi
 implicit none
 real, intent(in), dimension(ndim) :: xi
 
 !rhofunc = 2. + sin(2.*pi*xi(1)) 
 !rhofunc = max(2. - dot_product(xi,xi),1.)
 rhofunc = 1./exp(5.*((xi(1)-0.5)**2 + (xi(2)-0.5)**2))
 
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

logical function iallowed(xi,rhoi,np,massp,xp,hp,fac)
 use dimen_mhd, only:ndim,dndim
 implicit none
 integer, intent(in) :: np
 real, intent(in), dimension(ndim) :: xi
 real, intent(in) :: rhoi,massp,fac
 real, intent(in), dimension(ndim,np) :: xp
 real, intent(in), dimension(np) :: hp
 integer :: i
 real :: hi,hmin,r2

 hi = (massp/rhoi)**dndim
 !print*,'rhoi = ',rhoi,' hmin = ',hmin,' massp = ',massp
 iallowed = .true.
 !--exclude if any particle is closer than the allowed minimum separation
 overi: do i=1,np
    hmin = fac*hi !!(0.5*(hi + hp(i))) ! use average
    r2 = (xi(1)-xp(1,i))**2 + (xi(2)-xp(2,i))**2
    !r2 = dot_product(xi-xp(:,i),xi-xp(:,i))
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
