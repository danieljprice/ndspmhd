program kernelplot
 use kernel
 use options
 use loguns

 implicit none
 integer :: i,j,nkernels,nacross,ndown
 integer, dimension(10) :: iplotorder
 real :: q,q2
 real, dimension(0:ikern) :: dqkern
 logical :: samepage

 data iplotorder /0, 3, 9, 7, 4, 6, 4, 9, 0, 0/   ! order in which kernels are plotted
 nkernels = 6
 iprint = 6   ! make sure output from kernel setup goes to screen
 samepage = .false.
 nacross = 3
 ndown = 2

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'

 if (samepage) then
    call pgbegin(0,'?',1,1) 
    call pgenv(0.0,3.0,-3.5,1.7,0,1)
    call pglabel('r/h',' ','Smoothing kernels')
 else
    call pgbegin(0,'?',1,1) 
 endif

 do j=1,nkernels
    ikernel = iplotorder(j)
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
       call pgsch(0.6)
       !call pgenv(0.0,3.0,-3.5,1.7,1,1) !-3.5 1.7
       call danpgtile(j,nacross,ndown,0.0,3.2,-3.5,1.8,  &
                      'r/h',' ',TRIM(kernelname),1)
!!       call pgwnad(0.0,3.0,-3.5,1.7) ! pgwnad or pgswin
!!       call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)
!!       call pglabel('r/h',' ',TRIM(kernelname))
       call pgsci(1)
    endif

    call pgsls(1)
!
!--BEWARE! program will crash here if compiled incorrectly. This can happen
!  if compiled in double precision - e.g. if kernelND.f90 has been
!  compiled with the source code more recently than the last compilation for
!  this program.
!    
    print*,'   plotting kernel...(will crash here if double precision)'
    call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))
    call pgsls(2)
    call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern))
    call pgsls(3)
    call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
    call pgsci(1)
    call pgsls(1)
!
!--calculate dispersion relation for this kernel
!
    if (nkernels.eq.1) then
       call kernelstability1D
    endif

 enddo
 
 read*
!
!--plot stability separately
!
 if (nkernels.gt.1) then
    !!call pgsubp(nacross,ndown)   !--divide viewport into panels
    call pgsch(0.6)
    call pgpage
    do j=1,nkernels
       ikernel = iplotorder(j)
       call setkern
       dq2table = radkern2/real(ikern)
       do i=0,ikern
          q2 = i*dq2table
          q = sqrt(q2)
          dqkern(i) = q
       enddo
       call kernelstability1D
    enddo
 endif

 call pgend

end program kernelplot
