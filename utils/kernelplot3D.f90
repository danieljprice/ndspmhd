!--------------------------------------------------------------------
!
! Plotting program for kernels and kernel errors in 3D
!
!--------------------------------------------------------------------
program kernelplot3D
 use kernels
 use pagesetup, only:setpage2
 use legends,   only:legend
 use kernel_sums, only:get_kernel_sums,centre_on_particle,print_plot_options
 implicit none
 integer :: i,j,ikernel,ndim,ierr
 integer :: nkernels,nacross,ndown
 integer, dimension(10)   :: iplotorder
 real    :: hi,dh,hmin,hmax
 real    :: psep,rhoi,q2
 logical :: samepage,plotall,plot_cubic
 integer, parameter :: npts = 101 !241 !76 !51
 integer, parameter :: nsetups = 2
 real(kind=4), dimension(npts) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4), dimension(npts,nsetups) :: wnorm,gradwnorm,grad2wnorm
 integer, parameter :: nx = 32
 integer :: iplot
 real, parameter :: pi = 3.1415926536

 data iplotorder /0, 3, 20, 61, 62, 100, 94, 62, 63, 13/   ! order in which kernels are plotted
 !data iplotorder /0, 3, 24, 23, 42, 44, 64, 65, 14, 13/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 plotall = .false.

 !--options for plotit routine
 plot_cubic = .true.  ! plot cubic spline as a comparison
 centre_on_particle = .true.

 if (plotall) then
    nkernels = 99
 else
    nkernels = 6
 endif
 samepage = .false.
 nacross = 3
 ndown = 2
 ndim = 3

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'
 
 iplot = 1
 menu: do while (iplot /= 0)
 print*,'Select your plot (plots shown for all kernels and two different lattices): '
 call print_plot_options()
 print*,'Enter your selection now (0 = quit):'
 read*,iplot
 if (iplot <= 0) stop

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
 hmax = 2.8*psep
 dh   = (hmax-hmin)/real(npts - 1)
 
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
       call plot_kernel(j)
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

       call get_kernel_sums(iplot,nsetups,nx,ndim,npts,xplot,yplot,yplot2,&
                            wnorm,gradwnorm,grad2wnorm,psep,hmin,dh,radkern2)
       if (ikernel.eq.0) then
          cubic1 = yplot
          cubic2 = yplot2
       endif

       if (j.gt.0) then
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

subroutine plotit(j,iplot,xplot,yplot,yplot2,cubic1,cubic2)
 implicit none
 integer, intent(in) :: j,iplot
 real(kind=4), dimension(:), intent(in) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 
 select case(iplot)
 case(7,8)
    !ymin = -999.; ymax = 4999.
    ymin = -100.; ymax = 100.
    !ymin = -1.; ymax = 5.
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
    ymin = 0.99
    ymax = 1.01
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
    !call pgline(size(xplot),xplot,cubic1)
 endif
!--current kernel
 call pgsls(1)
 call pgline(size(xplot),xplot,yplot)

 call pgsci(3)
!--cubic spline
 if (plot_cubic) then
    call pgsls(4)
 !   call pgline(size(xplot),xplot,cubic2)
 endif
!--current kernel 
 call pgsls(2)
 call pgline(size(xplot),xplot,yplot2)

 call pgsci(1)
 call pgsls(1)
 
end subroutine plotit

subroutine plot_kernel(j)
 implicit none
 integer, intent(in) :: j
 real(kind=4) :: xmin,xmax,ymin,ymax,q
 integer :: i
 real(kind=4) :: dqkern(0:ikern),wstuff(0:ikern),wstuff2(0:ikern)
 character(len=20) :: ylabel

 xmin = 0. !minval(xplot)
 xmax = 3. !maxval(xplot)
 ymin = -0.8 !-2.2
 ymax = 0.8 !1.4
 ylabel = 'W, grad W, del^2 W'
!
!--setup x axis
!
 dq2table = radkern2/real(ikern)
 do i=0,ikern
    q2 = i*dq2table
    q = sqrt(q2)
    dqkern(i) = q
    !--this is h^ndim * W + h^(ndim+2) * del^2 W
    wstuff(i) = wij(i) + grgrwij(i) + (ndim - 1)/q*grwij(i)
    !--this is h^ndim * W + h^(ndim+2) * d^2 W/dh^2
    wstuff2(i) = 0.*wij(i) + ndim*(ndim+1)*wij(i) + 2.*(ndim+1)*q*grwij(i) + q**2*grgrwij(i)
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
 !--other term
 !call pgsci(4)
 !call pgsls(4)
 !call pgline(ikern+1,dqkern(0:ikern),wstuff(0:ikern))
 !call pgline(ikern+1,dqkern(0:ikern),wstuff2(0:ikern))
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
 ymin = 0.98
 ymax = 1.02
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

end program kernelplot3D
