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
 use dumpfiles, only:write_dump
!
!--define local variables
!            
 implicit none
 integer :: i,j,ipart,ix,iy,npartx,nparty,its,ixm1,iym1,incx,incy,ixmin,ixmax,iymin,iymax
 real :: totmass,rowmass,rowmassprev,rhofunc,rhomax,rhoprev,rhoi
 real :: dx,dy,massp,rij2,rijnew,rij,scale,dxrow
 real, dimension(ndim) :: xi,xj,xmid,dxpart
 integer, parameter :: maxits = 10000, nx = 100, ny = 100, nz = 100
 integer, external :: iget
 character(len=20) :: dumpfile
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

 print*,' enter npartx '
 read*,npartx
 nparty = npartx
 ntotal = npartx*nparty
 npart = ntotal
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

!
!--now setup particles on their initial cubic lattice
!
 dx = (xmax(1)-xmin(1))/real(npartx)
 dy = (xmax(2)-xmin(2))/real(nparty)

 ipart = 0
 do j=1,nparty
    xi(2) = (j-0.5)*dy
    do i=1,npartx
       xi(1) = (i-0.5)*dx
       ipart = ipart + 1
       x(:,ipart) = xi
    enddo
 enddo
 !volpart = product(xmax(:)-xmin(:))/real(ntotal)
 
! its = 0
! do while (its < maxits)
!    write(dumpfile,"(a,i6.6)") 'iteration',its
!    call write_dump(real(its),trim(dumpfile))!!
!
!    its = its + 1
!    print*,' ITERATION ',its
!    do i=1,npart
!       do j=i+1,npart
!          !--calculate separation of pair
!          dxpart(:) = x(:,i) - x(:,j)
!          rij2 = dot_product(dxpart,dxpart)
!          xmid(:) = 0.5*(x(:,i) + x(:,j))
!          rijnew = (massp/rhofunc(xmid))**dndim
!          if (abs(rij2 - rijnew**2)/rij2.lt.0.25) then
!            rij = sqrt(rij2)
 !            rijnew = max(rijnew,0.95*rij)
!             !--separate pair according to their density
!             x(:,i) = xmid(:) - 0.5*rijnew*dxpart/(rij+epsilon(rij2))
!             x(:,j) = xmid(:) + 0.5*rijnew*dxpart/(rij+epsilon(rij2))
!          endif
!       enddo
!    enddo
! enddo


!
!--now proceed with the iterations
!
 its = 0
 do while (its < maxits)
    if (mod(its,10).eq.0) then
       write(dumpfile,"(a,i6.6)") 'iteration',its
       call write_dump(real(its),trim(dumpfile))
    endif
    its = its + 1
    print*,' ITERATION ',its
    if (its .lt. maxits/4 .or. its.gt.3*maxits/4) then
       incx = 1
    else
       incx = -1
    endif
    if (its .lt. maxits/2) then
       incy = 1
    else
       incy = -1
    endif
    
    if (incx.eq.1) then
       ixmax = npartx
       ixmin = 1
    else
       ixmax = 1
       ixmin = npartx
    endif
    if (incy.eq.1) then
       iymax = nparty
       iymin = 1
    else
       iymax = 1
       iymin = nparty
    endif

    ipart = 0
    do iy=iymin,iymax,incy
       overx: do ix=ixmin,ixmax,incx
          if (incx.eq.1) then
             ixm1 = ix - 1
          else
             ixm1 = ix + 1
          endif
          ipart = ix + (iy-1)*npartx
          !--calculate current separation from previous particle in x
          xi(:) = x(:,ipart)
          if (incx.eq.1) then
             if (ix.gt.1) then
                xj(:) = x(:,iget(ixm1,iy,npartx,nparty))
             else
                xj(1) = xmin(1)
                xj(2) = xi(2)
             endif
          else
             if (ix.lt.npartx) then
                xj(:) = x(:,iget(ixm1,iy,npartx,nparty))
             else
                xj(1) = xmax(1)
                xj(2) = xi(2)
             endif          
          endif
          
          call adjust(xi,xj)
          !print*,' particle ',ipart,iget(ix-1,iy,npartx,nparty),' x = ',x(:,ipart),' xnew = ',xi(:)
          !read*
          
          if (xi(1).gt.xmax(1)) then
             print*,'end of row ',ix,xi(1)
             !exit overx
          endif

          if (incy.eq.1) then
             iym1 = iy - 1
          else
             iym1 = iy + 1
          endif          
          !--shift relative to previous particle in y direction
          if (incy.eq.1) then
             if (iy.gt.1) then
                xj(:) = x(:,iget(ix,iym1,npartx,nparty))
             else
                xj(1) = xi(1)
                xj(2) = xmin(2)
             endif
          else
             if (iy.lt.nparty) then
                xj(:) = x(:,iget(ix,iym1,npartx,nparty))
             else
                xj(1) = xi(1)
                xj(2) = xmax(2)
             endif          
          endif
          
          call adjust(xi,xj)
          !print*,' particle ',ipart,iget(ix,iy-1,npartx,nparty),' xnew = ',xi(:)
          !read*
          
          if (xi(2).gt.xmax(2)) then
             !print*,'end of column ',iy,xi(2)
             !exit overy
          endif
          
          x(:,ipart) = xi(:)

          !--at end of column, readjust the entire column
          if (incy.eq.1 .and. iy.eq.nparty) then
             !print*,'end of row ',ipart
             dxrow = (x(2,ipart) + 0.5*(massp/rhofunc(x(:,ipart)))**dndim) - xmin(2)
             scale = (xmax(2) - xmin(2))/dxrow
             !print*,'scale = ',scale
             !read*          
             do j = iymin,iymax,incy
                ipart = ix + (j-1)*npartx
                x(2,ipart) = xmin(2) + (x(2,ipart) - xmin(2))*scale
             enddo
          elseif (incy.eq.-1 .and. iy.eq.1) then
             dxrow = xmax(2) - (x(2,ipart) - 0.5*(massp/rhofunc(x(:,ipart)))**dndim)
             scale = (xmax(2) - xmin(2))/dxrow
             !print*,'scale = ',scale
             !read*
             do j = iymin,iymax,incy
                ipart = ix + (j-1)*npartx
                x(2,ipart) = xmax(2) - (xmax(2) - x(2,ipart))*scale
             enddo
          endif
       enddo overx
       !--at end of row, readjust the entire row
       !print*,'end of row ',ipart
       if (incx.eq.1) then
          dxrow = (x(1,ipart) + 0.5*(massp/rhofunc(x(:,ipart)))**dndim) - xmin(1)
          scale = (xmax(1) - xmin(1))/dxrow
          !print*,'scale = ',scale
          !read*
          do ix= ixmin,ixmax,incx
             ipart = ix + (iy-1)*npartx
             x(1,ipart) = xmin(1) + (x(1,ipart) - xmin(1))*scale
          enddo
       else
          dxrow = xmax(1) - (x(1,ipart) - 0.5*(massp/rhofunc(x(:,ipart)))**dndim)
          scale = (xmax(1) - xmin(1))/dxrow
          !print*,'scale = ',scale
          !read*
          do ix= ixmin,ixmax,incx
             ipart = ix + (iy-1)*npartx
             x(1,ipart) = xmax(1) - (xmax(1) - x(1,ipart))*scale
          enddo       
       endif

    enddo
 enddo
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
 
contains

subroutine adjust(xposi,xposj)
 implicit none
 real, intent(inout), dimension(ndim) :: xposi
 real, intent(in), dimension(ndim) :: xposj
 real, dimension(ndim) :: dxij,dr,xmid
 real :: rij,rijnew
 
 !--get separation vector from previous particle
 dxij(:) = xposi(:) - xposj(:)
 rij = sqrt(dot_product(dxij,dxij))
 dr(:) = dxij(:)/(rij + epsilon(rij))
 !print*,'dr = ',dxij(:),rij,dr(:)
 
 !--evaluate density at the midpoint between these two particles
 xmid(:) = 0.5*(xposi(:) + xposj(:))
 rijnew = (massp/rhofunc(xmid))**dndim
 !print*,'pos old = ',xposi(:),' pos new = ',xposj(:) + rijnew*dr(:),'rijnew = ',rijnew, ' rij = ',rij
 
 !--shift relative to previous x particle
 rijnew = min(rijnew,1.05*rij)
 rijnew = max(rijnew,0.95*rij)
 xposi(:) = xposj(:) + rijnew*dr(:)
 if (rijnew.lt.1.e-6) print*,' EEK ',rijnew

 return
end subroutine adjust

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

integer function iget(ix,iy,nx,ny)
 implicit none
 integer, intent(in) :: ix,iy,nx,ny
 
 iget = ix + (iy-1)*nx
 
end function

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
