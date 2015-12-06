!--------------------------------------------------------------------
!
! Module to compute various summations to do with SPH kernel errors
! in 3D. Used by kernelplot3D and kernelcalc3D
!
!--------------------------------------------------------------------
module kernel_sums
 implicit none
 public  :: get_kernel_sums, print_plot_options

 logical, public :: centre_on_particle = .true.

 private

contains

!-------------------------------------------
!+
!  print which sums are used for which plots
!+
!-------------------------------------------
subroutine print_plot_options()

 print*,' 1) kernels, derivatives and second derivatives as function of r/h'
 print*,' 2) kernel normalisation conditions as function of r/h '
 print*,' 3) kernel moments (R_xy) as function of r/h'
 print*,' 4) density as function of r/h'
 print*,' 5) gradw normalisation condition'
 print*,' 6) gradgradw normalisation condition'
 print*,' 7) del^2 rho/dh^2 as a function of r/h'
 print*,' 8) del^2 kernel stuff as a function of r/h'

end subroutine print_plot_options

!-------------------------------------------
!+
!  select which sum is used for which plot
!+
!-------------------------------------------
subroutine select_plot(iplot,sums,wsum,grad2sum,yploti)
 integer, intent(in)  :: iplot
 real,    intent(in)  :: sums(:),wsum,grad2sum
 real(4),    intent(out) :: yploti

 select case(iplot)
 case(8)
    yploti = sums(9)
 case(7)
    yploti = sums(7)
 case(6)
    yploti = grad2sum
 case(5)
    yploti = sums(4)
 case(4)
    yploti = wsum
 case(3)
    yploti = sums(2)
 case default
    yploti = 0.
 end select
 
end subroutine select_plot

!-------------------------------------------
!+
!  compute the various kernel summations
!+
!-------------------------------------------
subroutine get_kernel_sums(iplot,nsetups,nx,ndim,npts,xplot,yplot,yplot2, &
                           wnorm,gradwnorm,grad2wnorm,psep,hmin,dh,radkern2)
 use kernels, only:interpolate_kernels,interpolate_kernel_curl
 integer, intent(in)  :: iplot,nsetups,nx,ndim,npts
 real(kind=4), intent(out), dimension(npts) :: xplot,yplot,yplot2
 real(kind=4), intent(out), dimension(npts,nsetups) :: wnorm,gradwnorm,grad2wnorm
 real,     intent(in) :: psep,hmin,dh,radkern2
 real    :: q2,q,xi(3),rij2,mi,rhoi,rhoj
 real    :: hi,hi1,hi21,hnew
 real    :: dx(3),wsum,wabi,gradwi,gradgradwi,dum,rij,rij1,hfact
 real    :: volfrac,dzmax,dymax,erri,errmin
 real    :: dwdhi,d2wdh2i,gradhsum,omegai,func,dfdh,dhdrhoi,rhohi,gradwsum(3),grad2wsum
 real    :: xyzpart(3,2*nx**3+1)
 real    :: sums(10)
 integer :: i,n,nneigh,np,idim,imin,ipt
 integer :: ipart,maxiterations,iteration,isetup

 ipt = 0 
 over_lattices: do isetup=1,nsetups
!
!--set up uniform lattice of particles
! in box of 0->1
!
    call setpart(isetup,nx,np,ndim,xyzpart,ipart,dymax,dzmax)
    volfrac = dymax*dzmax
    if (ipart.ne.0) then
       maxiterations = 10
       xi(:) = xyzpart(:,ipart)
       !print*,' using particle ',ipart,' x,y,z = ',xi(1:ndim)
    else
       maxiterations = 10
       xi(:) = (nx/2 - 1)*psep + 0.25*psep
       !print*,' xi = ',xi
       !xi = 0.5
    endif
    maxiterations = 1
    
    !print*,' np  = ',np,' vol = ',volfrac
    mi   = 1./real(np)*volfrac
    errmin = huge(errmin)
    imin = 0
    rhoi = 1.

    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
       hfact = hi/psep
       if (abs(hfact-1.2).lt.0.01*dh) ipt = i
       !
       ! perform rho-h iterations to get actual h
       !
       its: do iteration=1,maxiterations
          !if (ipt.eq.i) hi = 1.2001545*psep

          hi1  = 1./hi
          hi21 = hi1*hi1
          !
          !--loop over neighbours, calculate sum
          !
          wsum = 0.
          gradhsum = 0.
          gradwsum(:) = 0.
          grad2wsum = 0.
          sums(:) = 0.
          nneigh = 0
          do n=1,np
             dx(:) = xi(:) - xyzpart(:,n)
             rij2 = dot_product(dx(1:ndim),dx(1:ndim))
             q2 = rij2*hi21
             if (q2.le.radkern2) then
                q      = sqrt(q2)
                rij    = sqrt(rij2)
                rij1   = 1./rij
                nneigh = nneigh + 1
                call interpolate_kernels(q2,wabi,dum,gradwi,gradgradwi)
                gradwi = 0.
                gradgradwi = 0.
                call interpolate_kernel_curl(q2,gradwi,gradgradwi)
                wabi   = wabi*hi1**ndim
                gradwi = gradwi*hi1**(ndim+1)
                gradgradwi = gradgradwi*hi1**(ndim+2)
                dwdhi  = -rij*gradwi*hi1 - ndim*wabi*hi1
                d2wdh2i = ndim*(ndim+1)*wabi*hi21 + 2.*(ndim+1)*q*gradwi*hi1 + gradgradwi*q2
                
                wsum = wsum + mi*wabi
                gradhsum = gradhsum + mi*dwdhi
                if (n.ne.ipart .and. q2.gt.0.) then
                   do idim=1,ndim
                      gradwsum(idim)  = gradwsum(idim) - mi/rhoi*dx(idim)*dx(idim)*rij1*gradwi
                   enddo
                   grad2wsum = grad2wsum + 0.5*mi/rhoi*dx(1)*dx(1)*gradgradwi
                   sums(1) = sums(1) + mi/rhoi*wabi*dx(1)*dx(1)/rij2
                   sums(2) = sums(2) + mi/rhoi*wabi*dx(1)*dx(2)/rij2
                   sums(3) = sums(3) + mi/rhoi*wabi*dx(1)*dx(3)/rij2
                   sums(4) = sums(4) - mi/rhoi*gradwi*dx(1)*dx(1)*rij1
                   sums(5) = sums(5) - mi/rhoi*gradwi*dx(1)*dx(2)*rij1
                   sums(6) = sums(6) - mi/rhoi*gradwi*dx(1)*dx(3)*rij1
                   rhoj = rhoi
                   sums(7) = sums(7) + rhoi*mi*(gradwi/rhoi**2 + gradwi/rhoj**2)*dx(1)*rij1
                   sums(8) = sums(8) + mi*d2wdh2i
                   sums(9) = sums(9) + mi*(gradgradwi + (ndim - 1)*gradwi*rij1) !hi1/sqrt(q2))
                   !if (iteration==maxiterations .and. isetup==2 .and. i==ipt) then
                   !print*,q2,gradgradwi/hi1**(ndim+2),mi*gradgradwi,mi*(ndim - 1)*gradwi*rij1
                   !endif
                endif
             endif
          enddo
          !
          !--Newton-Raphson stuff for density
          !
          rhoi     = wsum
          rhohi    = mi/(hi/hfact)**ndim
          dhdrhoi  = - hi/(ndim*wsum)
          omegai   = 1. - dhdrhoi*gradhsum
          gradhsum = 1./omegai
          func     = rhohi - wsum
          dfdh     = omegai/dhdrhoi
          
          erri = abs(wsum - 1.)
          
          if (maxiterations.gt.1) then
             hnew = hi - func/dfdh
             if (hnew < 0.8*hi) then
                hi = 0.8*hi
             elseif (hnew > 1.2*hi) then
                hi = 1.2*hi
             else
                hi = hnew
             endif
          endif
          sums(1:3) = ndim*sums(1:3)
          !sums(9) = sums(9)*psep**(ndim+2)/mi
         ! sums(9) = sums(9)*hi**(ndim+2)/mi

          wnorm(i,isetup) = wsum
          gradwnorm(i,isetup) = gradwsum(1)
          grad2wnorm(i,isetup) = 0.5*grad2wsum

          if (isetup.eq.1) then
             call select_plot(iplot,sums,wsum,0.5*grad2wsum,yplot(i))
          elseif (isetup.eq.2) then
             call select_plot(iplot,sums,wsum,0.5*grad2wsum,yplot2(i))
          endif
          if (i.eq.ipt .and. iteration.eq.maxiterations .and. isetup.eq.2) then
              !print*,'iteration ',iteration,' h/psep = ',hi/psep,' R = ',sqrt(radkern2)
             print*,'rho = ',wsum,' grad1 = ',sums(7),' yplot = ',yplot2(i),' del2rho = ',sums(9),hi1,' nneigh = ',nneigh
          endif
       enddo its

       if (erri.lt.errmin .and. xplot(i).gt.0.8 .and. xplot(i).lt.1.1) then
          imin = i
          errmin = erri
       endif
       !print*,i,' h = ',xplot(i),' wsum= ',wsum,' moments=',sums(1:3),' gradw = ',sums(4:6),&
       !         ' nneigh = ',nneigh,4./3.*pi*(radkern*hi)**3/mi
    enddo
    if (errmin.lt.0.1) print*,' best eta for lattice ',isetup,' at ',xplot(imin),' err = ',errmin
    
    !if (ipt.gt.1) write(isetup,*) abs(yplot(ipt)-1.)

 enddo over_lattices

end subroutine get_kernel_sums

!----------------------------------------------------------------------------
!+
!  routine to set up the particles on either a cubic or close-packed lattice
!+
!----------------------------------------------------------------------------
subroutine setpart(ilattice,nx,n,ndim,xyzpart,ipart,ymax,zmax)
 implicit none
 integer, intent(in) :: ilattice,nx,ndim
 integer, intent(out) :: n,ipart
 real, dimension(:,:), intent(out) :: xyzpart
 real, intent(out) :: ymax,zmax
 integer :: k,j,i,npartx,ny,nz,imin
 real :: xi,yi,zi,psep,xstart,ystart,zstart,r2,rmin
 real :: dx,dy,dz
 
 psep = 1./real(nx)
 
 n = 0
 ymax = 1.
 zmax = 1.
 xyzpart(:,:) = 0.
 if (ilattice.eq.1) then
 !--close-packed lattice
    dx = psep
    dy = 0.5*sqrt(3.)*psep
    dz = sqrt(6.)/3.*psep

    npartx = int(0.999/dx) + 1
    if (ndim.ge.2) then
       ny = int(0.999/dy) + 1
       !--adjust to exact multiples
       ny = 2*int(ny/2)
       ymax = ny*dy
    else
       ny = 1
    endif
    if (ndim.ge.3) then
       nz = int(0.999/dz) + 1
       nz = 3*int(nz/3)
       zmax = nz*dz
    else
       nz = 1
    endif
    
    do k=1,nz
       do j=1,ny
          ystart = dy/6.
          zstart = 0.5*dz
          xstart = 0.25*dx
          if (mod(k,3).eq.0) then  ! 3rd layer
             ystart = ystart + 2./3.*dy
             if (mod(j,2).eq.0) xstart = xstart + 0.5*dx
          elseif (mod(k,3).eq.2) then ! 2nd layer
             ystart = ystart + 1./3.*dy
             if (mod(j,2).eq.1) xstart = xstart + 0.5*dx
          elseif (mod(j,2).eq.0) then
             xstart = xstart + 0.5*dx
          endif
          do i = 1,npartx
             n = n + 1
             xyzpart(1,n) = (i-1)*dx + xstart
             if (ndim.ge.2) xyzpart(2,n) = (j-1)*dy + ystart
             if (ndim.eq.3) xyzpart(3,n) = (k-1)*dz + zstart
             !print*,n,' xyz = ',xyzpart(:,n)
          enddo
       enddo
    enddo
    !print*,' closepacked, setup ',n,' particles ',npartx,ny,nz,npartx*ny*nz

 else
 !--cubic lattice
    if (ndim.ge.3) then
       nz = nx
    else
       nz = 1
    endif
    if (ndim.ge.2) then
       ny = nx
    else
       ny = 1
    endif
    do k=1,nz
       zi = (k-1)*psep
       do j=1,ny
          yi = (j-1)*psep
          do i=1,nx
             xi = (i-1)*psep
       !print*,'x, y, z = ',xi,yi,zi
             n = n + 1
             xyzpart(1,n) = xi
             if (ndim.ge.2) xyzpart(2,n) = yi
             if (ndim.eq.3) xyzpart(3,n) = zi
          enddo
       enddo
    enddo
 endif
 
 rmin = huge(rmin)
 zi   = 0.
 yi   = 0.
 imin = 0
 do i=1,n
    xi = xyzpart(1,i) - 0.5
    if (ndim.ge.2) yi = xyzpart(2,i) - 0.5
    if (ndim.eq.3) zi = xyzpart(3,i) - 0.5
    r2 = xi*xi + yi*yi + zi*zi
    if (r2 .lt. rmin) then
       rmin = r2
       imin = i
    endif
 enddo
!
!--ipart is the particle location to be used in the density calculation
!  (if ipart=0 another position is chosen, not necessarily corresponding
!   to a particle in the setup)
!
 if (centre_on_particle) then
    !print*,' using particle ',imin,' at ',xyzpart(:,imin)
    ipart = imin
 else
    ipart = 0
 endif

end subroutine setpart

end module kernel_sums
