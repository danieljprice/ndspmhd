program kernelplot
 use kernel
 use options
 use loguns

 implicit none
 integer :: i, nkernels
 real :: q,q2
 real, dimension(0:ikern) :: dqkern
 logical :: samepage

 nkernels = 8
 iprint = 6   ! make sure output from kernel setup goes to screen
 samepage = .false.

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'

 call pgbegin(0,'?',1,1)
 if (samepage) then
    call pgenv(0.0,3.0,-3.5,1.7,0,1)
    call pglabel('r/h',' ','Smoothing kernels')
    call pgask(.false.)
 endif

 do ikernel=1,nkernels
    call setkern
!
!--setup x axis
!
    dq2table = radkern2/real(ikern)
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       dqkern(i) = q
    enddo
    if (samepage) then
       call pgsci(ikernel)
    else
       call pgenv(0.0,3.0,-3.5,1.7,0,1)
       call pglabel('r/h',' ',TRIM(kernelname))
       call pgask(.false.)
       call pgsci(1)
    endif

    call pgsls(1)
!!    print*,'plotting kernel...(will crash here if double precision)'
    call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))
    call pgsls(2)
    call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern))
    call pgsls(3)
    call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
    call pgsci(1)
    call pgsls(1)
    read*
 enddo

 call pgend

end program kernelplot
