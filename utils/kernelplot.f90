program kernelplot
 use debug
 use kernel
 use options
 use loguns
 use anticlumping
 use setup_params

 implicit none
 integer :: i,j,nkernels,nacross,ndown,ilinestyle
 integer, dimension(10) :: iplotorder
 real :: q,q2,xmin,xmax,ymin,ymax
 real, dimension(0:ikern) :: dqkern
 logical :: samepage
 character(len=50) :: string

 data iplotorder /0, 3, 10, 7, 5, 6, 4, 9, 0, 0/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 nkernels = 1
 iprint = 6   ! make sure output from kernel setup goes to screen
 idivBzero = 0
 eps = 0.4
 neps = 2
 hfact = 1.5
 samepage = .false.
 trace = .true.
 nacross = 1
 ndown = 1
 xmin = 0.0
 xmax = 3.2
 ymin = -3.5
 ymax = 1.7

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'

 if (samepage) then
    call pgbegin(0,'?',1,1) 
    !!call pgpap(11.7,0.6/sqrt(2.))  ! change paper size
    !!call pgenv(xmin,xmax,ymin,ymax,0,1)
    call pgsch(1.5)
    !!call pgenv(xmin,xmax,ymin,ymax,0,1)
    !!if (nkernels.gt.1) call pglabel('r/h',' ','Smoothing kernels')
    !!call pglabel('r/h',' ',' ')
    !!call pglabel('r/h','W, \(2266)W ',' ')
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
       call danpgtile(1,nacross,ndown,xmin,xmax,ymin,ymax,'r/h',' ',' ',0,1)
!!       if (nkernels.eq.1) call pglabel('r/h','W(r/h), \(2266)W(r/h) ',TRIM(kernelname))
    else
       call pgsch(0.6)
       !call pgenv(0.0,3.0,-3.5,1.7,1,1) !-3.5 1.7
       call danpgtile(j,nacross,ndown,xmin,xmax,ymin,ymax,  &
                      'r/h',' ',TRIM(kernelname),1,1)
!!       call pgwnad(0.0,3.0,-3.5,1.7) ! pgwnad or pgswin
!!       call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)
!!       call pglabel('r/h',' ',TRIM(kernelname))
       call pgsci(1)
    endif

!
!--BEWARE! program will crash here if compiled incorrectly. This can happen
!  if compiled in double precision - e.g. if kernelND.f90 has been
!  compiled with the source code more recently than the last compilation for
!  this program.
!    
    print*,'   plotting kernel...(will crash here if double precision)'
!
!--plot kernel
!
    call pgsls(1)
    call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))
!
!--kernel gradient
!    
    call pgsls(2)
    call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern))
!
!--second derivative
!
    call pgsls(3)
    call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
    call pgsci(1)
    call pgsls(1)

 enddo
 
 read*
!
!--plot stability separately
!
 print*,'-------------- plotting kernel stability ------------------'
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
       call kernelstability1D(j,nacross,ndown)
    enddo

 call pgend

end program kernelplot
!
!--draw a legend for different line styles
!  uses current line style and colour
!
subroutine legend(icall,text,hpos,vposin)
  implicit none
  integer, intent(in) :: icall
  character(len=*), intent(in) :: text
  real, intent(in) :: vposin
  real, dimension(2) :: xline,yline
  real :: xch, ych, xmin, xmax, ymin, ymax
  real :: vspace, hpos, vpos

  call pgstbg(0)           ! opaque text to overwrite previous
!
!--set horizontal and vertical position and spacing
!  in units of the character height
!
  vspace = 1.5  ! (in units of character heights)
!  hpos = 0.4   ! distance from edge, in fraction of viewport
  vpos = vposin + (icall-1)*vspace  ! distance from top, in units of char height

  call pgqwin(xmin,xmax,ymin,ymax) ! query xmax, ymax
  call pgqcs(4,xch,ych) ! query character height in x and y units 

  yline = ymax - ((vpos - 0.5)*ych)
  xline(1) = xmin + hpos*(xmax-xmin)
  xline(2) = xline(1) + 3.*xch

  call pgline(2,xline,yline)            ! draw line segment
!!--make up line style if > 5 calls (must match actual line drawn)
  if (icall.eq.1) then
     call pgpt(2,xline,yline,17) !mod(icall,5)+1)
     call pgpt(1,0.5*(xline(1)+xline(2)),yline(1),17) !mod(icall,5)+1)
  endif  
  call PGTEXT(xline(2) + 0.5*xch,yline(1)-0.25*ych,trim(text))

!!  call pgmtxt('T',vpos,0.1,0.0,trim(text))  ! write text

  call pgstbg(-1) ! reset text background to transparent

end subroutine legend
