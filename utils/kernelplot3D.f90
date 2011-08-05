program kernelplot3D
 use kernels
 use pagesetup, only:setpage2
 use legends,   only:legend
 implicit none
 integer :: i,j,n,ikernel,ndim,ierr,imin
 integer :: nkernels,nacross,ndown,nneigh
 integer :: ipart,maxiterations,iteration,isetup
 integer, dimension(10)   :: iplotorder
 real    :: q2,xi,yi,zi,rij2,mi,rhoi
 real    :: hi,hi1,hi21,hi31,hi41,dh,hmin,hmax
 real    :: psep,dx,dy,dz,wsum,wabi,gradwi,dum,rij,rij1,hfact
 real    :: volfrac,dzmax,dymax,erri,errmin
 real    :: dwdhi,gradhsum,omegai,func,dfdh,dhdrhoi,rhohi
 logical :: samepage,plotall,plot_moments, plot_kernels
!! character(len=50) :: text
 integer, parameter :: npts = 101 !241 !76 !51
 real(kind=4), dimension(npts) :: xplot,yplot,yplot2,cubic1,cubic2
 integer, parameter :: nx = 32
 integer :: np,ipt
 real, dimension(3,2*nx**3+1) :: xyzpart
 real, parameter :: pi = 3.1415926536
 real, dimension(10) :: sums

 !data iplotorder /0, 3, 51, 52, 53, 59, 64, 65, 14, 13/   ! order in which kernels are plotted
 data iplotorder /0, 3, 42, 43, 42, 44, 64, 65, 14, 13/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 plotall = .false.
 plot_moments = .true. ! either plot kernel moments or density/normalisation conditions
 plot_kernels = .false.

 if (plotall) then
    nkernels = 70
 else
    nkernels = 1
 endif
 samepage = .false.
 nacross = 2
 ndown = 2
 ndim = 3

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'

 if (samepage) then
    call pgbegin(0,'?',1,1) 
    if (nacross.eq.2 .and. ndown.eq.1) call pgpap(11.7,0.6/sqrt(2.))  ! change paper size
    call pgsch(1.)
 else
    call pgpap(6.,1.)
    !call pgpap(11.5,0.75)
    call pgbegin(0,'?',1,1)
 endif
 call pgslw(1)

 psep = 1./real(nx)
 rhoi = 1.
 hmin = 1.0*psep
 hmax = 1.8*psep
 dh   = (hmax-hmin)/real(npts - 1)
 ipt = 0
 
 if (plot_kernels) then
    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
    enddo
    do j=1,nkernels
       if (j.gt.0 .and. j.le.size(iplotorder)) ikernel = iplotorder(j)
       call setkern(ikernel,ndim,ierr)
       print*,'plotting ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2
       call plot_kernel(j,xplot)
    enddo
 endif
 
 do j=0,nkernels
    if (plotall) then
       ikernel = j
    else
       if (j.gt.0) ikernel = iplotorder(j)
    endif
    call setkern(ikernel,ndim,ierr)
    print "(60('-'))"
    print*,'Using kernel ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2

    over_lattices: do isetup=1,2
!
!--set up uniform lattice of particles
! in box of 0->1
!
    call setpart(isetup,nx,np,xyzpart,ipart,dymax,dzmax)
    volfrac = dymax*dzmax
    if (ipart.ne.0) then
       maxiterations = 10
       xi = xyzpart(1,ipart)
       yi = xyzpart(2,ipart)
       zi = xyzpart(3,ipart)
       !print*,' using particle ',ipart,' x,y,z = ',xi,yi,zi
    else
       maxiterations = 1
       xi = (nx/2 - 1)*psep + 0.25*psep
       print*,' xi = ',xi
       !xi = 0.5
       yi = xi
       zi = xi
    endif
    !print*,' np  = ',np,' vol = ',volfrac
    mi   = 1./real(np)*volfrac
    errmin = huge(errmin)
    imin = 0

    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
       if (abs(xplot(i)-1.2).lt.0.01*dh) ipt = i
       hfact = hi/psep
       
       its: do iteration=1,maxiterations
          hi1  = 1./hi
          hi21 = hi1*hi1
          hi31 = hi21*hi1
          hi41 = hi21*hi21
          !
          !--loop over neighbours, calculate sum
          !
          wsum = 0.
          gradhsum = 0.
          sums(:) = 0.
          nneigh = 0
          do n=1,np
             dx = xi - xyzpart(1,n)
             dy = yi - xyzpart(2,n)
             dz = zi - xyzpart(3,n)
             rij2 = dx*dx + dy*dy + dz*dz
             q2 = rij2*hi21
             if (q2.le.radkern2) then
                rij    = sqrt(rij2)
                rij1   = 1./rij
                nneigh = nneigh + 1
                call interpolate_kernels(q2,wabi,gradwi,dum,dum)
                wabi   = wabi*hi31
                gradwi = gradwi*hi41
                dwdhi  = -rij*gradwi*hi1 - ndim*wabi*hi1
                
                wsum = wsum + mi*wabi
                gradhsum = gradhsum + mi*dwdhi
                if (n.ne.ipart .and. q2.gt.0.) then
                   sums(1) = sums(1) + mi/rhoi*wabi*dx*dx/rij2
                   sums(2) = sums(2) + mi/rhoi*wabi*dx*dy/rij2
                   sums(3) = sums(3) + mi/rhoi*wabi*dx*dz/rij2
                endif
             endif
          enddo
          !
          !--Newton-Raphson stuff for density
          !
          rhohi    = mi/(hi/hfact)**ndim
          dhdrhoi  = - hi/(ndim*wsum)
          omegai   = 1. - dhdrhoi*gradhsum
          gradhsum = 1./omegai
          func     = rhohi - wsum
          dfdh     = omegai/dhdrhoi
          
          erri = abs(wsum - 1.)
          
          hi = hi - func/dfdh
          sums(:) = ndim*sums(:)
          if (isetup.eq.1) then
             yplot(i) = wsum 
             !yplot(i) = log10(erri)
             if (ikernel.eq.0) cubic1(i) = yplot(i)
             !yplot(i) = sums(2)
          endif
          if (isetup.eq.2) then
             yplot2(i) = wsum
             !yplot2(i) = log10(erri)
             !yplot2(i) = sums(2)
             if (ikernel.eq.0) cubic2(i) = yplot2(i)
          endif
          if (ikernel.eq.65 .and. i.eq.ipt) print*,'hfact = ',hfact,': iteration ',iteration,' hi = ',hi,' rho = ',5.*wsum     
          !yplot2(i) = sums(2)
       enddo its
       
       if (erri.lt.errmin .and. xplot(i).gt.0.8 .and. xplot(i).lt.1.1) then
          imin = i
          errmin = erri
       endif
       print*,i,' h = ',xplot(i),' wsum= ',wsum,' moments=',sums(1:3),' nneigh = ',nneigh,4./3.*pi*(radkern*hi)**3/mi
    enddo
    if (errmin.lt.0.1) print*,' best eta for '//trim(kernelname)//' (lattice ',isetup,') at ',xplot(imin),' err = ',errmin
    
    if (ipt.gt.1) write(isetup,*) ikernel,abs(yplot(ipt)-1.)

    enddo over_lattices

    if (j.gt.0) then
       if (ipt.gt.0) then
         if (abs(yplot(ipt)-1.).lt.abs(cubic1(ipt)-1.)) then
            print*,' BEATS cubic on cubic by ',abs(cubic1(ipt)-1.)/abs(yplot(ipt)-1.),' at ',xplot(ipt)
         endif
         if (abs(yplot2(ipt)-1.).lt.abs(cubic2(ipt)-1.)) then
            print*,' BEATS cubic on closep by ',abs(cubic2(ipt)-1.)/abs(yplot2(ipt)-1.),' wsum = ',wsum
         endif
         print*,' wsum (cubic,closep) = ',5.*yplot(ipt),5.*yplot2(ipt),' for cubic spline = ',5.*cubic1(ipt),5.*cubic2(ipt)
       endif
       if (plot_kernels) then
          call plotit(nkernels+j,xplot,yplot,yplot2,cubic1,cubic2)
       else
          call plotit(j,xplot,yplot,yplot2,cubic1,cubic2)
       endif
       !read*
    endif
 enddo
 call pgend
 read*

contains

subroutine plot_kernel(j,xplot)
 implicit none
 integer, intent(in) :: j
 real(kind=4), dimension(:), intent(in) :: xplot
 real(kind=4) :: xmin,xmax,ymin,ymax,q
 integer :: i
 real(kind=4), dimension(0:ikern) :: dqkern
 character(len=20) :: ylabel

 xmin = 0. !minval(xplot)
 xmax = 3. !maxval(xplot)
 ymin = -2.2
 ymax = 1.4
 ylabel = 'W(norm)'
!
!--setup x axis
!
 dq2table = radkern2/real(ikern)
 do i=0,ikern
    q2 = i*dq2table
    q = sqrt(q2)
    dqkern(i) = q
 enddo
 print*,'xmin,max = ',xmin,xmax
 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'r/h',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
 call pgmtxt('t',-2.0_4,0.93_4,1.0_4,trim(kernelname))
 !--kernel function
 call pgsci(1)
 call pgsls(1)
 call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern)) 
 !--gradient
 call pgsls(2)
 call pgsci(2)
 call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern)) 
 !--2nd deriv
 call pgsci(3)
 call pgsls(3)
 call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
 !--restore settings
 call pgsci(1)
 call pgsls(1)

end subroutine plot_kernel

subroutine plotit(j,xplot,yplot,yplot2,cubic1,cubic2)
 implicit none
 integer, intent(in) :: j
 real(kind=4), dimension(:), intent(in) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = 1. !minval(xplot)
 xmax = 1.4 !maxval(xplot)
 !ymin = -0.08
 !ymax = 0.12
 !ymin = minval(yplot2) - 0.005
 !ymax = maxval(yplot2) + 0.005
 !if (j.eq.1) then
 ymin = 0.997
 ymax = 1.001
 !endif
 ylabel = 'W(norm)'
 ylabel = 'R_{xy}'

 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'h/r',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
 call pgmtxt('t',-2.0_4,0.96_4,1.0_4,trim(kernelname))
!--cubic spline
 call pgsci(2)
! call pgsls(4)
! call pgline(size(xplot),xplot,cubic1)
!--current kernel 
 call pgsls(1)
 call pgline(size(xplot),xplot,yplot)

 call pgsci(2)
!--cubic spline
! call pgsls(4)
! call pgline(size(xplot),xplot,cubic2)
!--current kernel 
 call pgsls(2)
 call pgline(size(xplot),xplot,yplot2)


 call pgsci(1)
 call pgsls(1)
 
end subroutine plotit

subroutine setpart(ilattice,nx,n,xyzpart,ipart,ymax,zmax)
 implicit none
 integer, intent(in) :: ilattice,nx
 integer, intent(out) :: n,ipart
 real, dimension(:,:), intent(out) :: xyzpart
 real, intent(out) :: ymax,zmax
 integer :: k,j,i,npartx,ny,nz,imin
 real :: xi,yi,zi,psep,xstart,ystart,zstart,r2,rmin
 
 psep = 1./real(nx)
 
 n = 0
 ymax = 1.
 zmax = 1.
 if (ilattice.eq.2) then
 !--close-packed lattice
    dx = psep
    dy = 0.5*sqrt(3.)*psep
    dz = sqrt(6.)/3.*psep

    npartx = int(0.999/dx) + 1
    ny = int(0.999/dy) + 1
    nz = int(0.999/dz) + 1
    
    !--adjust to exact multiples
    ny = 2*int(ny/2)
    nz = 3*int(nz/3)
    ymax = ny*dy
    zmax = nz*dz
    
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
             xyzpart(2,n) = (j-1)*dy + ystart
             xyzpart(3,n) = (k-1)*dz + zstart
             !print*,n,' xyz = ',xyzpart(:,n)
          enddo
       enddo
    enddo
   ! print*,' closepacked, setup ',n,' particles ',npartx,ny,nz,npartx*ny*nz

 else
 !--cubic lattice
    do k=1,nx
       zi = (k-1)*psep
       do j=1,nx
          yi = (j-1)*psep
          do i=1,nx
             xi = (i-1)*psep
       !print*,'x, y, z = ',xi,yi,zi
             n = n + 1
             xyzpart(1,n) = xi
             xyzpart(2,n) = yi
             xyzpart(3,n) = zi
          enddo
       enddo
    enddo
 endif
 
 rmin = huge(rmin)
 do i=1,n
    xi = xyzpart(1,i) - 0.5
    yi = xyzpart(2,i) - 0.5
    zi = xyzpart(3,i) - 0.5
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
 if (plot_moments) then
    ipart = 0
 else
    print*,' particle ',imin,' at ',xyzpart(:,imin)
    ipart = imin
 endif

end subroutine setpart

end program kernelplot3D
