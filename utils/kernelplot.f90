program kernelplot
 use kernels
 implicit none
 integer :: i,j,ikernel,ndim
 integer :: nkernels,nacross,ndown,nepszero,ipos,i1,i2
 integer, dimension(10) :: iplotorder
 real :: q,q2,xmin,xmax,ymin,ymax,epszero,hpos,vpos,rhoi
 real, dimension(0:ikern) :: dqkern
 logical :: samepage
!! character(len=50) :: text

 data iplotorder /1, 6, 31, 21, 16, 16, 16, 16, 16, 16/   ! order in which kernels are plotted
 !!iplotorder = 0 ! override data statement if all the same kernel
 nkernels = 4
 ianticlump = 0
 epszero = 0.0
 nepszero = 4
 samepage = .false.
 nacross = 4
 ndown = 1
 xmin = 0.0
 xmax = 3.2
!! ymin = 0.01
 ymin = -3.5  !!3.5
!! ymax = 2.4
 ymax = 2.7
 ipos = 1
 ndim = 1

 print*,'welcome to kernel city, where the grass is green and the kernels are pretty...'
 print*,'plotting kernels normalised in ',ndim,' dimensions'

 if (samepage) then
    call pgbegin(0,'?',1,1) 
    if (nacross.eq.2 .and. ndown.eq.1) call pgpap(11.7,0.6/sqrt(2.))  ! change paper size
    !!call pgenv(xmin,xmax,ymin,ymax,0,1)
    call pgsch(1.)
    !!call pgenv(xmin,xmax,ymin,ymax,0,1)
    !!if (nkernels.gt.1) call pglabel('r/h',' ','Smoothing kernels')
    !!call pglabel('r/h',' ',' ')
    !!call pglabel('r/h','W, \(2266)W ',' ')
 else
    call pgbegin(0,'?',1,1) 
    !!call pgpap(5.85,2./sqrt(2.))
 endif
 call pgslw(2)

 eps = epszero
 neps = nepszero

 do j=1,nkernels
    if (ianticlump.ne.0) then !! .and. j.lt.7) then
       eps = eps + 0.4
       print*,'eps = ',eps, ' neps = ',neps
!    else
!       if (j.eq.7) eps = 0.2
!       eps = eps + 0.2
    endif
    
    ikernel = iplotorder(j)
    call setkern(ikernel,ndim)
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
       hpos = 0.7
       vpos = 3.0
!       if (j.lt.7) then
       call danpgtile(1,nacross,ndown,xmin,xmax,ymin,ymax,'r/h',' ',' ',0,1)
       !
       !--plot legend
       !
!        call legend(1,'W',hpos,vpos)
!        call pgsls(2)
!        call legend(2,'\(2266)W',hpos,vpos)
!       call pgsls(1)
!	call legend(0,'\ge = 0.0-1.0',hpos-0.125,vpos+4.5)      
!        write(text,"(a,i1)") 'n = ',neps
!        call legend(0,text,hpos-0.125,vpos+6.0)

!       else
!        call danpgtile(2,nacross,ndown,xmin,xmax,ymin,ymax,'r/h',' ',' ',0,1)
!       !
!       !--plot legend
!       !
!        call legend(1,'W',hpos,vpos)
!        call pgsls(2)
!        call legend(2,'\(2266)W',hpos,vpos)
!        call pgsls(1)       
!        write(text,"(a,f3.1)") '\ge = ',eps
!        call legend(0,text,hpos-0.125,vpos+4.5)
!	call legend(0,'n = 3-5',hpos-0.125,vpos+6.0)

!       endif
!!       if (nkernels.eq.1) call pglabel('r/h','W(r/h), \(2266)W(r/h) ',TRIM(kernelname))
       
    else
       !!call pgsch(1.5)
       if (mod(j,nacross*ndown).eq.1 .or. nacross*ndown.eq.1) call pgpage
       !call pgenv(0.0,3.0,-3.5,1.7,1,1) !-3.5 1.7
       !!ymax = maxval(wij(1:ikern))*1.5
       if (ikernel.eq.101) ymax = maxval(abs(wij(1:ikern)))*1.5
       if (ikernel.eq.102) ymax = maxval(wij(1:ikern)/dqkern(1:ikern)**2)*1.5
       print*,'max = ',ymax
       if (ikernel.eq.0) kernelname = 'density'
       call danpgtile(j,nacross,ndown,xmin,xmax,ymin,ymax,  &
                      'r/h',' ',TRIM(kernelname),0,1)
!!       call pgwnad(0.0,3.0,-3.5,1.7) ! pgwnad or pgswin
!!       call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)
!!       call pglabel('r/h',' ',TRIM(kernelname))
       call pgsci(1)
       ipos = 1
    endif
!
!--test kernel normalisation if 1D
!  assume h=1, h=1.2*psep, rho=1, equispaced particles
!  so 2h = 2.4, neighbours at 1/1.2 and 2.0/1.2
!
   if (ndim.eq.1) then
      i1 = int((1./1.2)**2*ddq2table)
      i2 = int((2./1.2)**2*ddq2table)
      rhoi = (wij(0) + 2.*(wij(i1) + wij(i2)))/1.2
      print*,'rhoi = ',rhoi
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
    if (ikernel.eq.101) wij(0:ikern) = -wij(0:ikern)
    if (ikernel.eq.102) wij(1:ikern) = wij(1:ikern)/dqkern(1:ikern)**2
    call pgline(ikern+1,dqkern(0:ikern),wij(0:ikern))
    if (ianticlump.ne.0) then
!       !!call pgsls(j+2)
       call pgline(ikern+1,dqkern(0:ikern),wijaniso(0:ikern))
    endif
    !!call legend(ipos,trim(kernelname),0.5,3.0)
    if (ikernel.eq.101) then
       call pgsls(2)
       call pgline(ikern,dqkern(1:ikern),1./dqkern(1:ikern))
       call legend(ipos,'1/r',0.78,3.0)
    elseif (ikernel.eq.102) then
       call pgsls(2)
       call pgline(ikern+1,dqkern(0:ikern),1./dqkern(0:ikern)**2)
       call legend(ipos,'1/r\u2',0.78,3.0)
    endif
!
!--kernel gradient
!    
    if (ikernel.lt.100) then
       call pgsls(2)
       call pgline(ikern+1,dqkern(0:ikern),grwij(0:ikern))
       !!call pgline(ikern+1,dqkern(0:ikern),dqkern(0:ikern)*grwij(0:ikern))
       if (ianticlump.ne.0) then
   !       call pgsls(j+2)
          !!write(text,"(a,f3.1,a,i1)") '\ge = ',eps,', n = ',neps
          !!call legend(ipos,text,0.5,3.0)
          call pgline(ikern+1,dqkern(0:ikern),grwijaniso(0:ikern))
         !! call pgline(ikern+1,dqkern(0:ikern),dqkern(0:ikern)*grwijaniso(0:ikern))    
       endif
    endif
!
!--second derivative
!
    if (ikernel.lt.100) then
       call pgsls(3)
       call pgline(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
       if (ianticlump.ne.0) then
          call pgline(ikern+1,dqkern(0:ikern),grgrwijaniso(0:ikern))    
       endif
       call pgsls(1)
!    call pgsci(1)
    endif
 enddo
 
!! call pgend
!
!--plot stability separately
!
 print*,'-------------- plotting kernel stability ------------------'
    !!call pgsch(0.6)
    call pgpage
!!    if (nkernels.eq.1) call pgpap(5.85,1./sqrt(2.))  ! change paper size

    eps = epszero
    neps = nepszero
    do j=1,nkernels
       if (ianticlump.ne.0) then
!          if (j.eq.2) then
!	     eps = epszero
!	     call pgsls(2)
!	  endif
          eps = eps + 0.4
       endif
       call setkern(iplotorder(j),ndim)
       call kernelstability1D(j,nacross,ndown,eps,neps)
    enddo

 call pgend

end program kernelplot
!
!--draw a legend for different line styles
!  uses current line style and colour
!
subroutine legend(icall,text,hpos,vposin)
  implicit none
  integer, intent(inout) :: icall
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

!!--make up line style if > 5 calls (must match actual line drawn)
!   if (icall.eq.3) then
!     call pgpt(2,xline,yline,17) !mod(icall,5)+1)
!     call pgpt(1,0.5*(xline(1)+xline(2)),yline(1),17) !mod(icall,5)+1)
!     !!call pgline(2,xline,yline)            ! draw line segment
   if (icall.gt.0) then
     call pgline(2,xline,yline)            ! draw line segment
   endif
   icall = icall + 1
!
!--write text next to line segment
!  
  call PGTEXT(xline(2) + 0.5*xch,yline(1)-0.25*ych,trim(text))

!!  call pgmtxt('T',vpos,0.1,0.0,trim(text))  ! write text

  call pgstbg(-1) ! reset text background to transparent

end subroutine legend
