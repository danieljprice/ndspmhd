program kernelplot3D
 use kernels
 use pagesetup, only:setpage2
 use legends,   only:legend
 implicit none
 integer :: i,j,n,ikernel,ndim,ierr,imin
 integer :: nkernels,nacross,ndown,nneigh
 integer :: ipart,maxiterations,iteration,isetup
 integer, dimension(10)   :: iplotorder
 real    :: q2,xi(3),rij2,mi,rhoi,rhoj
 real    :: hi,hi1,hi21,dh,hmin,hmax,hnew
 real    :: psep,dx(3),wsum,wabi,gradwi,gradgradwi,dum,rij,rij1,hfact
 real    :: volfrac,dzmax,dymax,erri,errmin
 real    :: dwdhi,gradhsum,omegai,func,dfdh,dhdrhoi,rhohi,gradwsum(3),grad2wsum
 logical :: samepage,plotall,plot_cubic
 logical :: centre_on_particle
!! character(len=50) :: text
 integer, parameter :: npts = 101 !241 !76 !51
 integer, parameter :: nsetups = 2
 real(kind=4), dimension(npts) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4), dimension(npts,nsetups) :: wnorm,gradwnorm,grad2wnorm
 integer, parameter :: nx = 32
 integer :: np,ipt,idim,iplot
 real, dimension(3,2*nx**3+1) :: xyzpart
 real, parameter :: pi = 3.1415926536
 real, dimension(10) :: sums

 data iplotorder /0, 2, 3, 61, 62, 63, 94, 62, 63, 13/   ! order in which kernels are plotted
 !data iplotorder /0, 3, 24, 23, 42, 44, 64, 65, 14, 13/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 plotall = .false.
 centre_on_particle = .true.

 !--options for plotit routine
 plot_cubic = .true.  ! plot cubic spline as a comparison

 if (plotall) then
    nkernels = 99
 else
    nkernels = 6
 endif
 samepage = .false.
 nacross = 3
 ndown = 2
 ndim = 1

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'
 
 iplot = 1
 menu: do while (iplot /= 0)
 print*,'Select your plot (plots shown for all kernels and two different lattices): '
 print*,' 1) kernels, derivatives and second derivatives as function of r/h'
 print*,' 2) kernel normalisation conditions as function of r/h '
 print*,' 3) kernel moments (R_xy) as function of r/h'
 print*,' 4) density as function of r/h'
 print*,' 5) gradw normalisation condition'
 print*,' 6) gradgradw normalisation condition'
 print*,' 7) del^2 rho as a function of r/h'
 print*,'Enter your selection now (-ve to NOT plot kernels as well; 0 = quit) '
 read*,iplot

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
 hmin = 0.8*psep
 hmax = 1.8*psep
 dh   = (hmax-hmin)/real(npts - 1)
 ipt = 0
 
 if (iplot==1) then
    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
    enddo
    do j=1,nkernels
       if (j.gt.0 .and. j.le.size(iplotorder)) ikernel = iplotorder(j)
       if (plotall) ikernel = j
       call setkern(ikernel,ndim,ierr)
       print*,'plotting ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2
       call plot_kernel(j,xplot)
    enddo
 else
 
 do j=0,nkernels
    if (plotall) then
       ikernel = j
    else
       if (j.gt.0) ikernel = iplotorder(j)
    endif
    call setkernels(ikernel,ikernel,ndim,ierr,ierr)
!    call setkern(ikernel,ndim,ierr)
    print "(60('-'))"
    print*,'Using kernel ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2

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
                rij    = sqrt(rij2)
                rij1   = 1./rij
                nneigh = nneigh + 1
                call interpolate_kernels(q2,wabi,dum,gradwi,gradgradwi)
                wabi   = wabi*hi1**ndim
                gradwi = gradwi*hi1**(ndim+1)
                gradgradwi = gradgradwi*hi1**(ndim+2)
                dwdhi  = -rij*gradwi*hi1 - ndim*wabi*hi1
                
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
                   sums(8) = sums(8) + mi*gradgradwi
                endif
             endif
          enddo
          !
          !--Newton-Raphson stuff for density
          !
          if (i.eq.ipt .and. iteration.eq.maxiterations) then
              print*,'iteration ',iteration,' h/psep = ',hi/psep,' R = ',sqrt(radkern2)
              print*,'rho = ',wsum,' grad1 = ',sums(7),' del2rho = ',sums(8)
          endif

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

          wnorm(i,isetup) = wsum
          gradwnorm(i,isetup) = gradwsum(1)
          grad2wnorm(i,isetup) = 0.5*grad2wsum

          if (isetup.eq.1) then
             call select_plot(iplot,ikernel,sums,wsum,0.5*grad2wsum,yplot(i),cubic1(i))
          elseif (isetup.eq.2) then
             call select_plot(iplot,ikernel,sums,wsum,0.5*grad2wsum,yplot2(i),cubic2(i))
          endif
       enddo its

       if (ikernel.eq.13) then
           print*,'eta =',(hmin + (i-1)*dh)/psep,'h = ',hi/psep,' GOT sums(8) = ',sums(8)
       endif
       
       if (erri.lt.errmin .and. xplot(i).gt.0.8 .and. xplot(i).lt.1.1) then
          imin = i
          errmin = erri
       endif
       !print*,i,' h = ',xplot(i),' wsum= ',wsum,' moments=',sums(1:3),' gradw = ',sums(4:6),&
       !         ' nneigh = ',nneigh,4./3.*pi*(radkern*hi)**3/mi
    enddo
    if (errmin.lt.0.1) print*,' best eta for '//trim(kernelname)//' (lattice ',isetup,') at ',xplot(imin),' err = ',errmin
    
    if (ipt.gt.1) write(isetup,*) ikernel,abs(yplot(ipt)-1.)

    enddo over_lattices

    if (j.gt.0) then
       if (.false. .and. ipt.gt.0) then
         if (abs(yplot(ipt)-1.).lt.abs(cubic1(ipt)-1.)) then
            print*,' BEATS cubic on cubic by ',abs(cubic1(ipt)-1.)/abs(yplot(ipt)-1.),' at ',xplot(ipt)
         endif
         if (abs(yplot2(ipt)-1.).lt.abs(cubic2(ipt)-1.)) then
            print*,' BEATS cubic on closep by ',abs(cubic2(ipt)-1.)/abs(yplot2(ipt)-1.),' wsum = ',wsum
         endif
         print*,' wsum (cubic,closep) = ',5.*yplot(ipt),5.*yplot2(ipt),' for cubic spline = ',5.*cubic1(ipt),5.*cubic2(ipt)
       endif
       select case(iplot)
       case(3)
          call plot_moments(j,xplot,yplot,yplot2,cubic1,cubic2)
       case(2)
          call plot_normalisations(nkernels+j,xplot,wnorm,gradwnorm,grad2wnorm)
          !call plotit(nkernels+j,xplot,yplot,yplot2,cubic1,cubic2)
       case default
          call plotit(j,iplot,xplot,yplot,yplot2,cubic1,cubic2)
       end select
       !read*
    endif
 enddo
 
 endif
 call pgend
 
 enddo menu

contains

subroutine select_plot(iplot,ikernel,sums,wsum,grad2sum,yploti,cubici)
 integer, intent(in)  :: iplot,ikernel
 real,    intent(in)  :: sums(:),wsum,grad2sum
 real,    intent(out) :: yploti,cubici

 select case(iplot)
 case(7)
    yploti = sums(8)
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
 if (ikernel.eq.0) cubici = yploti
 
end subroutine select_plot
 

subroutine plotit(j,iplot,xplot,yplot,yplot2,cubic1,cubic2)
 implicit none
 integer, intent(in) :: j,iplot
 real(kind=4), dimension(:), intent(in) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 
 select case(iplot)
 case(7)
    ymin = -999.0
    ymax = 2999.
    ylabel = 'del^2 \rho'
 case(6)
    ymin = 0.9
    ymax = 1.1
    ylabel = '\nabla^2 W_{norm}'
 case(5)
    ymin = 0.9
    ymax = 1.1
    ylabel = '\nabla W_{norm}'
 case default
    ymin = 0.95
    ymax = 1.08
    ylabel = 'W(norm)'
 end select

 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'\eta',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
 call pgmtxt('t',-2.0_4,0.96_4,1.0_4,trim(kernelname))
!--cubic spline
 call pgsci(2)
 if (plot_cubic) then
    call pgsls(4)
    call pgline(size(xplot),xplot,cubic1)
 endif
!--current kernel
 call pgsls(1)
 call pgline(size(xplot),xplot,yplot)

 call pgsci(3)
!--cubic spline
 if (plot_cubic) then
    call pgsls(4)
    call pgline(size(xplot),xplot,cubic2)
 endif
!--current kernel 
 call pgsls(2)
 call pgline(size(xplot),xplot,yplot2)

 call pgsci(1)
 call pgsls(1)
 
end subroutine plotit

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
 ylabel = 'W, grad W, del^2 W'
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

subroutine plot_normalisations(j,xplot,wnormi,grwnorm,gr2wnorm)
 implicit none
 integer, intent(in) :: j
 real(kind=4), dimension(:), intent(in) :: xplot
 real(kind=4), dimension(:,:), intent(in) :: wnormi,grwnorm,gr2wnorm
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 ymin = 0.95
 ymax = 1.05
 ylabel = 'W_{norm}'

 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'h/\gDx',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)

 call pgmtxt('t',-2.0_4,0.96_4,1.0_4,trim(kernelname))
 do i=1,size(wnormi(1,:))
    call pgsci(i+1)

   !--normalisation of W (solid line)
    call pgsls(1)
    call pgline(size(xplot),xplot,wnormi(:,i))

   !--normalisation of grad W (dashed line)
    call pgsls(2)
    call pgline(size(xplot),xplot,grwnorm(:,i))

   !--normalisation of del^2 W (dotted line)
    call pgsls(4)
    call pgline(size(xplot),xplot,gr2wnorm(:,i))
    print*,i,' gr2wnorm = ',gr2wnorm(1:10,i)
 enddo

 call pgsci(1)
 call pgsls(1)
 
end subroutine plot_normalisations

subroutine plot_moments(j,xplot,yplot,yplot2,cubic1,cubic2)
 implicit none
 integer, intent(in) :: j
 real(kind=4), dimension(:), intent(in) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 ymin = -0.05
 ymax = 0.12
 ylabel = 'R_{xy}' 

 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'h/\gDx',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
 call pgmtxt('t',-2.0_4,0.96_4,1.0_4,trim(kernelname))
!--cubic spline
 call pgsci(2)
 if (plot_cubic) then
    call pgsls(4)
    call pgline(size(xplot),xplot,cubic1)
 endif
!--current kernel
 call pgsls(1)
 call pgline(size(xplot),xplot,yplot)

 call pgsci(3)
!--cubic spline
 if (plot_cubic) then
    call pgsls(4)
    call pgline(size(xplot),xplot,cubic2)
 endif
!--current kernel 
 call pgsls(2)
 call pgline(size(xplot),xplot,yplot2)

 call pgsci(1)
 call pgsls(1)
 
end subroutine plot_moments

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

end program kernelplot3D
