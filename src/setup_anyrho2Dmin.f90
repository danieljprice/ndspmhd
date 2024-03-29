module functions
 use dimen_mhd, only:ndim
 use part, only:x
 integer :: iset
 real :: massp
 
contains

real function distfunc(xi)
 implicit none
 real, intent(in), dimension(ndim) :: xi
 real :: rhoi,rcut,sum,rhofunc,r2
 integer :: idim,j

 rhoi = rhofunc(xi)
 rcut = (massp/rhoi)**(1./ndim)

 sum = 0.
 do j=1,iset-1
    r2 = 0.
    do idim=1,ndim
       r2 = r2 + (xi(idim) - x(idim,j))**2
    enddo
    sum = sum + rcut**4/r2 + r2
 enddo
 distfunc = sum !+ rcut*rcut
 
end function distfunc 

subroutine ddistfunc(xi,grad)
 implicit none
 real, intent(in), dimension(ndim) :: xi
 real, intent(out), dimension(ndim) :: grad
 real :: rhoi,rcut,rhofunc,r2
 integer :: j,idim
 
 rhoi = rhofunc(xi)
 rcut = (massp/rhoi)**(1./ndim)

 grad = 0.
 do j=1,iset-1
    r2 = 0.
    do idim=1,ndim
       r2 = r2 + (xi(idim) - x(idim,j))**2
       grad(idim) = grad(idim) + 2.*(xi(idim) - x(idim,j))
    enddo

    do idim=1,ndim
       grad(idim) = grad(idim)*(1. - rcut**4/(r2*r2))
    enddo
 enddo
 
end subroutine ddistfunc

end module functions
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
 use random,    only:ran1,ran2
 use functions, only:distfunc,ddistfunc,iset,massp
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
 real :: totmass,rowmass,rowmassprev,rhofunc,fac,drho
 real :: myran,rhomax,xi(ndim),dx,dy,rhoi,rhoprev,dxtry,dxmin
 real :: dxstep,xmincell,ymincell,fret,t1,t2
 real, dimension(ndim) :: grad
 logical, external :: iallowed
 logical, parameter :: userandom = .true.,adaptmesh=.false.,usemin=.false.
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
 
 call cpu_time(t1)
 
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
 call sobseq(-2,xi)
 if (userandom .and. .not.adaptmesh) then
    fac = 0.85*sqrt(1./0.834)
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
    if (usemin) then
       iset = 1
       x(1,iset) = 0.5 !xmin(1)
       x(2,iset) = 0.5 !xmin(2)
       do ipart = 2,ntotal
          iset = ipart
         ! print*,'here ',xi,rhofunc(xi)
          !call drhofunc(xi,grad)
          !print*,'grad = ',grad
          xi(1) = xmin(1) + ran2(iseed) !xi(1) + dxmin*(ran2(iseed)-0.5)
          xi(2) = xmin(2) + ran2(iseed) !xi(2) + dxmin*(ran2(iseed)-0.5)
          call dfpmin(xi,ndim,1.e-6,itry,fret,distfunc,ddistfunc)
          print*,ipart,' x = ',xi,' fret = ',fret,' its = ',itry
          print*,'func = ',distfunc(xi(:)),rhofunc(xi),(massp/rhofunc(xi))**(1./ndim)
          call ddistfunc(xi,grad)
          print*,'deriv = ',grad(:)
          x(:,iset) = xi(:)
       enddo
       ipart = iset
    else
 
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
          if (userandom) then
             do itry=1,2*nxpts
             
                call sobseq(ndim,xi)
                xi(1) = xmincell + dx*xi(1)
                xi(2) = ymincell + dy*xi(2)
                rhoi = rhofunc(xi)
                if (iallowed(xi,rhoi,ipart,massp,x(:,1:ipart),hh(1:ipart),fac)) then
                   ipart = ipart + 1
                   x(:,ipart) = xi
                   print*,ipart,xi,rhoi
                   hh(ipart) = (massp/rhoi)**dndim
                endif
             enddo 
       
          else
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
          endif
!          icell = nextcell(icell)
!       enddo
!       print*,' ncells this level = ',ncellslevel,' npart= ',ipart
!       print*,' npart = ',ipart
    enddo
    
    endif
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
       !if (mod(itry,10000).eq.0) then
       !   fac = fac*0.9999
       !   print*,fac
       !   if (fac.lt.0.85) ntotal = ipart
       !else
       !   fac = 0.9*sqrt(1./0.834)
       !endif
       if (itry.gt.1000000) ntotal = ipart
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
 call cpu_time(t2)
 print*,'npart =',npart, ' total time for setup = ',t2-t1,'s' 
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

end subroutine setup

real function rhofunc(xi)
 use dimen_mhd, only:ndim
 use setup_params, only:pi
 implicit none
 real, intent(in), dimension(ndim) :: xi
 
 !rhofunc = dot_product(xi-0.25,xi-0.25)
 !rhofunc = 2.
 rhofunc = 2. + sin(2.*pi*xi(1)) 
 !rhofunc = max(2. - dot_product(xi,xi),1.)
 !rhofunc = 1./exp(5.*((xi(1)-0.5)**2 + (xi(2)-0.5)**2))
 
end function rhofunc

subroutine drhofunc(xi,grad)
 use dimen_mhd, only:ndim
 implicit none
 real, intent(in), dimension(ndim) :: xi
 real, intent(out), dimension(ndim) :: grad
 integer :: i
 
 do i=1,ndim
    grad(i) = 2.*(xi(i) - 0.25)
 enddo
 
end subroutine drhofunc

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
