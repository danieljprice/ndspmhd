program kernelplot3D
 use kernels
 use pagesetup, only:setpage2
 use legends,   only:legend
 implicit none
 integer :: i,j,k,n,ikernel,ndim,ierr
 integer :: nkernels,nacross,ndown,nneigh
 integer, dimension(10)   :: iplotorder
 real    :: q2,xmin,xmax,ymin,ymax,xi,yi,zi,rij2,mi,rhoi,hi,hi1,hi21,hi31,dh,hmin,hmax
 real    :: psep,dx,dy,dz,wsum,wabi,gradwi,dum
 logical :: samepage
!! character(len=50) :: text
 integer, parameter :: npts = 50
 real, dimension(npts) :: xplot,yplot
 integer, parameter :: nx = 32
 integer, parameter :: np = nx*nx*nx
 real, dimension(3,np+1) :: xyzpart

 data iplotorder /0, 2, 3, 10, 64, 65, 14, 13, 16, 16/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 nkernels = 4
 samepage = .false.
 nacross = 4
 ndown = 1
 ndim = 3

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'
 
 if (samepage) then
    call pgbegin(0,'?',1,1) 
    if (nacross.eq.2 .and. ndown.eq.1) call pgpap(11.7,0.6/sqrt(2.))  ! change paper size
    call pgsch(1.)
 else
    call pgpap(8.5,0.5)
    call pgbegin(0,'?',1,1)
 endif
 call pgslw(1)

!
!--set up cubic lattice of particles
! in box of 0->1
!
 psep = 1./real(nx-1) 
 n  = 0
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
 
 xi = 0.5
 yi = 0.5
 zi = 0.5
 xyzpart(1,np+1) = xi
 xyzpart(2,np+1) = yi
 xyzpart(3,np+1) = zi
 mi   = 1./real(np)
 rhoi = 1.
 hmin = psep
 hmax = 2.*psep
 dh   = (hmax-hmin)/real(npts - 1)
 
 do j=1,nkernels
    ikernel = iplotorder(j)
    call setkern(ikernel,ndim,ierr)
    print*,'Using '//trim(kernelname)//' kernel, r^2 = ',radkern2,ndim

    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
       hi1  = 1./hi
       hi21 = hi1*hi1
       hi31 = hi21*hi1
       !
       !--loop over neighbours, calculate sum
       !
       wsum = 0.
       nneigh = 0
       do n=1,np
          dx = xi - xyzpart(1,n)
          dy = yi - xyzpart(2,n)
          dz = zi - xyzpart(3,n)
          rij2 = dx*dx + dy*dy + dz*dz
          q2 = rij2*hi21
          if (q2.le.radkern2) then
             nneigh = nneigh + 1
             call interpolate_kernels(q2,wabi,gradwi,dum,dum)
             wsum = wsum + mi*wabi*hi31
          endif
       enddo
       yplot(i) = wsum
       print*,' h = ',xplot(i),' wsum = ',yplot(i),' nneigh = ',nneigh
    enddo

    call plotit(j,xplot,yplot)
    read*
 enddo


contains

subroutine plotit(j,xplot,yplot)
 implicit none
 integer, intent(in) :: j
 real, dimension(:), intent(in) :: xplot,yplot

 xmin = minval(xplot)
 xmax = maxval(xplot)
 ymin = 0.8 !minval(yplot)
 ymax = 1.2 !maxval(yplot)
 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'r/h','W(norm)', &
               trim(kernelname),0,1,&
               0.,0.,0.,0.,0.,0.,.false.,.true.)
 call pgmtxt('t',-2.0,0.93,1.0,trim(kernelname))
 call pgline(size(xplot),xplot,yplot)
 
end subroutine plotit

end program kernelplot3D
