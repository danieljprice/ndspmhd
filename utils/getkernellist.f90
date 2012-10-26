program getkernellist
 use kernels, only:setkern,kernelname
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
       write(lu,"(a)") trim(kernelname)
    else
       write(lu,"(a,i3)") 'nokernel ',i
    endif
 enddo
 close(lu)

end program getkernellist
