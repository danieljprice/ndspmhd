!--------------------------------------------------------------------
!
!  Like kernelplot3D but does not do any plotting. Instead, we
!  simply return errors for use in kernel breeding
!
!--------------------------------------------------------------------
program kernelcalc3D
 use kernels,     only:kernelname,setkern,radkern2,verbose,setkernels,radkern,write_kernel_table
 use kernel_sums, only:get_kernel_sums,print_plot_options
 integer, parameter :: npts = 101
 integer, parameter :: nsetups = 2
 integer, parameter :: maxkernels = 121
 real, dimension(npts) :: xplot,yplot,yplot2
 real, dimension(npts,nsetups) :: wnorm,gradwnorm,grad2wnorm
 integer, parameter :: nx = 32
 integer :: iplot,nargs,ndim,i1,i2
 character(len=32) :: string
 character(len=len(kernelname)) :: kname0
 real :: L2err(nsetups),Linf(nsetups),exact,erri
 real, parameter :: pi = 3.1415926536

 nargs = command_argument_count()
 ndim = 3
 if (nargs < 2) then
    print "(/,a)",'Usage: kernelcalc ikernel iplot'
    print "(/,a)",' Where iplot is one of the following:'
    call print_plot_options()
    print "(/,a)",' press any key for kernel list:'
    read*
    ! list all available kernels
    kname0 = ''
    verbose = .false.
    write_kernel_table = .true.
    do i=0,maxkernels
       call setkern(i,ndim,ierr)
       if (ierr == 0 .and. (trim(kernelname) /= trim(kname0))) then
          print "(i3,': ',a)",i,trim(kernelname)
       endif
       if (i==0) kname0 = kernelname
    enddo
    stop
 endif
 ! extract command line arguments for ikernel and iplot
 call get_command_argument(1,string)
 read(string,*) ikernel
 if (ikernel < 0) stop 'need ikernel >= 0'
 call get_command_argument(2,string)
 read(string,*) iplot
 if (iplot <= 0) stop 'need iplot > 0'

 print*,'computing kernels normalised in ',ndim,' dimensions'

 psep = 1./real(nx)
 rhoi = 1.
 hmin = 1.0*psep
 hmax = 2.*psep
 dh   = (hmax-hmin)/real(npts - 1)
 
 write_kernel_table = .false.
 call setkernels(ikernel,ikernel,ndim,ierr,ierr)
 print*,'Using kernel ',ikernel,': '//trim(kernelname)//' kernel, r^2 = ',radkern2

 print "(60('-'))"

 call get_kernel_sums(iplot,nsetups,nx,ndim,npts,xplot,yplot,yplot2,&
                      wnorm,gradwnorm,grad2wnorm,psep,hmin,dh,radkern2)
 
 select case(iplot)
 case (7,8)
    exact = 0.
 case default
    exact = 1.
 end select

 i1 = 10
 i2 = 11
 open(i1,file='errors-cp.out',status='replace')
 open(i2,file='errors-cubic.out',status='replace')
 do isetup=1,nsetups
    L2err(isetup) = 0.
    Linf(isetup)  = 0.
    do i=1,npts
       xneighi = 4./3.*pi*(xplot(i)*radkern)**3
       if (isetup==1) then
          erri = abs(yplot(i)-exact)*(xneighi/32.)
          write(i1,*) xplot(i),erri,xneighi,yplot(i)
       else
          erri = abs(yplot2(i)-exact)*(xneighi/32.)   
          write(i2,*) xplot(i),erri,xneighi,yplot2(i)
       endif
       L2err(isetup) = L2err(isetup) + erri**2
       Linf(isetup)  = max(Linf(isetup),erri)
    enddo
    L2err(isetup) = sqrt(L2err(isetup)/npts)
 enddo
 close(i1)
 close(i2)
 print*, ' max error at h = ',xplot(maxloc(yplot)),' and h = ',xplot(maxloc(yplot2))
 print*,' L2 Error (close packed) = ',L2err(1),' Linf = ',Linf(1)
 print*,' L2 error (cubic)        = ',L2err(2),' Linf = ',Linf(2)
 
 print*,' score:',sqrt(L2err(1)**2 + L2err(2)**2)

end program kernelcalc3D
