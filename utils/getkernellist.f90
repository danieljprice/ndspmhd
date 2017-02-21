program getkernellist
 use kernels, only:setkern,kernelname,radkern2
 implicit none
 integer :: i,ndim,ierr
 integer, parameter :: lu = 1
 
 ndim = 1
 open(unit=lu,file='kernels.list',status='replace')
 do i=1,120
    print*,i
    kernelname=''
    call setkern(i,ndim,ierr)
    if (len_trim(kernelname).gt.0) then
       if (radkern2 <= 4. .and. .not. (kernelname=='M_4 cubic spline')) then
       write(lu,"(i3,':',a)") i,trim(kernelname)
       endif
    else
       write(lu,"(a,i3)") 'nokernel ',i
    endif
 enddo
 close(lu)

end program getkernellist
