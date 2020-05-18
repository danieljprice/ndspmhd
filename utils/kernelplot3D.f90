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
 logical :: samepage,plotall,plot_cubic,phantom_paper,use_doublehump_limits
 integer, parameter :: npts = 101 !241 !76 !51
 integer, parameter :: nsetups = 1
 real(kind=4), dimension(npts) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4), dimension(npts,nsetups) :: wnorm,gradwnorm,grad2wnorm,grad2wnormb
 integer, parameter :: nx = 32
 integer :: iplot
 real, parameter :: pi = 3.1415926536

 !data iplotorder /0, 105, 59, 104, 103, 3, 94, 62, 63, 13/   ! order in which kernels are plotted
 !data iplotorder /0, 3, 106, 107, 94, 108, 64, 65, 14, 13/   ! order in which kernels are plotted
 !data iplotorder /0, 2, 3, 61, 62, 63, 64, 65, 14, 13/   ! order in which kernels are plotted
 !data iplotorder /41, 42, 43, 70, 71, 72, 64, 65, 14, 13/   ! order in which kernels are plotted
 !data iplotorder /4,93,5,95,60,3,92,43,2,20/   ! order in which kernels are plotted
 data iplotorder /20,94,99,0,3,200,31,90,3,100/   ! order in which kernels are plotted
 !data iplotorder /100,11,98,42,09,64,61,24,30,33/
 !!iplotorder = 0 ! override data statement if all the same kernel
 plotall = .false.

 !--options for plotit routine
 plot_cubic = .false.  ! plot cubic spline as a comparison
 centre_on_particle = .true.
 phantom_paper = .false.
 use_doublehump_limits = .false.

 if (plotall) then
    nkernels = 99
 else
    nkernels = 6
 endif
 samepage = .false.
 nacross = 3
 ndown = 2
 ndim = 3

 print*,'Welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'Plotting kernels normalised in ',ndim,' dimensions'
 
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
    if (phantom_paper) call pgpap(10.,1.0)
    call pgsch(1.1)
 endif
 call pgslw(1)

 psep = 1./real(nx)
 rhoi = 1.
 hmin = 1.0*psep
 hmax = 2.8*psep
 dh   = (hmax-hmin)/real(npts - 1)
 
 if (iplot==1 .or. iplot==10) then
    do i=1,npts
       hi = hmin + (i-1)*dh
       xplot(i) = hi/psep
    enddo
    do j=1,nkernels
       if (j.gt.0 .and. j.le.size(iplotorder)) ikernel = iplotorder(j)
       if (plotall) ikernel = j
       call setkern(ikernel,ndim,ierr,r=2.)
       print*,'plotting ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2
       if (iplot==10) then
          call plot_softening(j)
       else
          call plot_kernel(j)
       endif
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
       print*,'Using kernel ',ikernel,': '//trim(kernelname)//' kernel, R = ',sqrt(radkern2)

       call get_kernel_sums(iplot,nsetups,nx,ndim,npts,xplot,yplot,yplot2,&
                            wnorm,gradwnorm,grad2wnorm,grad2wnormb,psep,hmin,dh,radkern2)
       if (ikernel.eq.0) then
          cubic1 = yplot
          cubic2 = yplot2
       endif

       if (j.gt.0) then
          select case(iplot)
          case(3)
             call plot_moments(j,xplot,yplot,yplot2,cubic1,cubic2)
          case(2)
             call plot_normalisations(nkernels+j,xplot,wnorm,gradwnorm,grad2wnorm,grad2wnormb)
             !call plotit(nkernels+j,xplot,yplot,yplot2,cubic1,cubic2)
          case default
             print*,yplot(:)
             call plotit(j,iplot,nsetups,xplot,yplot,yplot2,cubic1,cubic2)
          end select
          !read*
       endif
    enddo
 endif
 call pgend
 
 enddo menu

contains

!----------------------------------------------------
!+
!  General plotting routine for a variety of things
!+
!----------------------------------------------------
subroutine plotit(j,iplot,nsetups,xplot,yplot,yplot2,cubic1,cubic2)
 implicit none
 integer, intent(in) :: j,iplot,nsetups
 real(kind=4), dimension(:), intent(in) :: xplot,yplot,yplot2,cubic1,cubic2
 real(kind=4) :: xmin,xmax,ymin,ymax
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 
 select case(iplot)
 case(11)
    ymin = 0.9
    ymax = 1.1
    ylabel = 'graddivv'
 case(8,9)
    ymin = -4999.; ymax = 4999.
    !ymin = -1000.; ymax = 1000.
    !ymin = -1.; ymax = 5.
    ylabel = 'del^2 \rho'
 case(7)
    ymin = -0.001
    ymax = 0.001
    ylabel = '0.5 \Sigma x_{ab} y_{ab} \nabla^2 W_{ab}'
 case(6)
    ymin = 0.7
    ymax = 1.1
    ylabel = '\nabla^2 W_{norm}'
 case(5)
    ymin = 0.98
    ymax = 1.02
    ylabel = '\nabla W_{norm}'
 case default
    ymin = 0.99
    ymax = 1.01
    ylabel = 'W(norm)'
    !ymin = 0.999
    !ymax = 1.001
    !ylabel = 'W(norm)'
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
 if (nsetups >= 2) then
!--current kernel 
    call pgsls(2)
    call pgline(size(xplot),xplot,yplot2)
 endif

 call pgsci(1)
 call pgsls(1)
 
end subroutine plotit

!---------------------------------------
!+
!  Plot kernel functions
!+
!---------------------------------------
subroutine plot_kernel(j)
 implicit none
 integer, intent(in) :: j
 real(kind=4) :: xmin,xmax,ymin,ymax,q
 integer :: i
 real(kind=4) :: dqkern(0:ikern),wstuff(0:ikern),wstuff2(0:ikern)
 character(len=20) :: ylabel

 xmin = 0. !minval(xplot)
 xmax = 2.99 !maxval(xplot)
 if (use_doublehump_limits) then
    ymin = 0.
    ymax = 0.39 
 else
    ymin = -0.8 ! -5.19 !-0.8 !-2.2
    ymax = 1.4 !5.19 !0.8 !1.4
 endif
 !ylabel = 'W, grad W, del^2 W'
 ylabel = 'f(q)'
!
!--setup x axis
!
 dq2table = radkern2/real(ikern)
 do i=0,ikern
    q2 = i*dq2table
    q = sqrt(q2)
    dqkern(i) = q
    !--this is h^ndim * W + h^(ndim+2) * del^2 W
    !wstuff(i) = wij(i) + grgrwij(i) + (ndim - 1)/q*grwij(i)
    !--this is h^ndim * W + h^(ndim+2) * d^2 W/dh^2
    !wstuff2(i) = 0.*wij(i) + ndim*(ndim+1)*wij(i) + 2.*(ndim+1)*q*grwij(i) + q**2*grgrwij(i)
    wstuff(i) = -2.*grwij(i)/q
    wstuff2(i) = grgrwij(i) + (ndim-1)*grwij(i)/q
 enddo
 call setpage2(j,nacross,ndown,xmin,xmax,ymin,ymax,'r/h',trim(ylabel), &
               trim(kernelname),0,1,&
               0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
 call pgmtxt('t',-2.0_4,0.93_4,1.0_4,trim(kernelname))
 
 call pgsci(1)
 call pgsls(1)
 !--kernel function
 call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))
 
 if (use_doublehump_limits) return
 !--gradient
 call pgsls(2)
 !call pgsci(2)
 call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern)) 
 !--2nd deriv
 !call pgsci(3)
 call pgsls(4)
 call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
 !--other term
 call pgsci(4)
 call pgsls(4)
 call pgline(ikern+1,dqkern(0:ikern),wstuff(0:ikern))
 call pgsci(5)
 call pgsls(5)
 call pgline(ikern+1,dqkern(0:ikern),wstuff2(0:ikern))
 !--restore settings
 call pgsci(1)
 call pgsls(1)

end subroutine plot_kernel

!------------------------------------
!+
!  Plot softening kernel functions
!+
!------------------------------------
subroutine plot_softening(j)
 use legends, only:legend
 integer, intent(in) :: j
 real(kind=4) :: xmin,xmax,ymin,ymax(3),q
 integer :: i,irow,icall,ix,iy,iref
 real(kind=4) :: dqkern(0:ikern),fref(0:ikern),phiref(0:ikern)
 character(len=20) :: ylabel(3)

 xmin = 0.
 xmax = 2.99
 ymin = 0.001
 ymax = (/1.09,1.96,1.4/)
 ylabel = (/'w(q)    ','\phi(q) ','\phi''(q)'/)
!
!--setup x axis
!
 dq2table = radkern2/real(ikern)
 do i=0,ikern
    q2 = i*dq2table
    q = sqrt(q2)
    dqkern(i) = q
    phiref(i) = 1./q
    fref(i) = 1./q2
 enddo
 do irow=1,3
    iy = (j-1)/nacross + 1
    ix = j - ((j-1)/nacross)*nacross
    iref = ix + (irow-1)*3 + (iy-1)*nacross*3

    call setpage2(iref,nacross,ndown*3,xmin,xmax,ymin,ymax(irow),'r/h',trim(ylabel(irow)), &
                  trim(kernelname),0,1,&
                  0._4,0._4,0._4,0._4,0._4,0._4,.false.,.true.)
    if (irow==1) call pgmtxt('t',-1.5_4,0.93_4,1.0_4,trim(kernelname))

    call pgsci(1)
    call pgsls(1)
    
    icall = 1
    select case(irow)
    case(1)
       !--kernel function
       call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))  
    case(3)
       !--force
       call pgline(ikern+1,dqkern(0:ikern),fsoft(0:ikern))
       !--1/r2
       call pgsls(4)
       call pgline(ikern+1,dqkern(0:ikern),fref(0:ikern))
       if (ix==3 .and. iy==1) call legend(icall,'1/r^2',0.5_4,2.4_4)
    case(2)
       !--potential
       call pgline(ikern+1,dqkern(0:ikern),-potensoft(0:ikern)) 
       !--1/r
       call pgsls(4)
       call pgline(ikern+1,dqkern(0:ikern),phiref(0:ikern))
       if (ix==3 .and. iy==1) call legend(icall,'1/r',0.5_4,2.4_4)
    end select
    call pgsls(1)
 enddo

end subroutine plot_softening

!---------------------------------------
!+
!  Plot kernel normalisation functions
!+
!---------------------------------------
subroutine plot_normalisations(j,xplot,wnormi,grwnorm,gr2wnorm,gr2wnormb)
 implicit none
 integer, intent(in) :: j
 real(kind=4), dimension(:), intent(in) :: xplot
 real(kind=4), dimension(:,:), intent(in) :: wnormi,grwnorm,gr2wnorm,gr2wnormb
 real(kind=4) :: xmin,xmax,ylim(2)
 character(len=20) :: ylabel

 xmin = minval(xplot)
 xmax = maxval(xplot)
 ylim = (/0.98,1.02/)
 !ylim = (/0.9,1.1/)
 ylabel = 'W_{norm}'

 call setpage2(j,nacross,ndown,xmin,xmax,ylim(1),ylim(2),'h/\gDx',trim(ylabel), &
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
    print*,i,' gr2wnorm = ',gr2wnorm(:,i)

   !--normalisation of del^2 W using Brookshaw method (dash-dotted line)
    !call pgsls(5)
    !call pgline(size(xplot),xplot,gr2wnormb(:,i))
    !print*,i,' gr2wnormb = ',gr2wnormb(:,i)
 enddo

 call pgsci(1)
 call pgsls(1)
 
end subroutine plot_normalisations

!---------------------------------------
!+
!  Plot kernel moments
!+
!---------------------------------------
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
