subroutine kernelstability1D(iplot,nacrossin,ndownin,eps,neps)
  use pagesetup, only:setpage2
  use kernels, only:kernelname
  use legends, only:legend
  implicit none
  integer, parameter :: ny = 100, nkx = 100, npart = 40, ncont = 40
  integer, intent(in) :: iplot, nacrossin, ndownin, neps
  real, intent(in) :: eps
  integer :: i,j,ipart,ikx,ipos
  integer :: iplotpos, iloop, nplots, nacross, ndown
  real, parameter :: pi = 3.1415926536
  real, dimension(nkx,ny) :: dat
  real, dimension(ny) :: yaxis,tterm1,tterm2,tterm3,tterm4,rhosum
  real, dimension(nkx) :: kxarray
  real, dimension(npart) :: x
  real, dimension(ncont) :: levels
  real, dimension(6) :: trans
  real :: xmin, psep, pmass, rhozero, przero, gamma, cs2, R
  real :: ymin, ymax, dy, dkx, kxmin, kxmax, h, kx
  real :: datmin, datmax, dcont, omegasq, omegasq1D
  real :: charheight, hpos, vpos, term1,term2,term3,term4,rhosumi
  character(len=5) :: string,labely
  character(len=50) :: text,title
  logical :: negstress, contours
  common /terms/ term1,term2,term3,term4,rhosumi
  
  cs2 = 1.0
  rhozero = 1.0
  przero = 1.0
  gamma = 1.0  ! isothermal
  psep = 1.0
  pmass = rhozero*psep
!
!--set various options
!  
  contours = .false.     ! plot whole dispersion relation or just kx=0
  negstress = .false.    ! plot vs h or R
  R = 1.0       ! R=1 gives usual hydrodynamics, R < 0 gives negative stress
  h = 1.2*psep   ! value of smoothing length
  ikx = 1  !!nkx/2        ! frequency for cs vs R or h plots
  nplots = 1
  nacross = nacrossin
  ndown = ndownin
!  nacross = int(nplots/2) + 1
!  ndown = nplots/nacross
!  nacross = 1
!  ndown = 1

!
!--iplotpos is plot position on page (= iplot if kernelstability being called
!  multiple times with different kernels, but different if multiple plots are
!  being drawn on the page within this subroutine).
!  
  iplotpos = iplot
  title = trim(kernelname)
!
!--set up an array of particle positions to sum over
!
  xmin = -npart/2*psep
  do i=1,npart
     x(i) = xmin + i*psep
  enddo
  ipart = npart/2
!
!--set up an array of smoothing length values for the y axis
!
  if (negstress) then   ! y axis is negative stress parameter R
     ymin = -10.
     ymax = 1.
     labely = 'R'
  else                  ! y axis is h
!     ymin = 0.5*psep
!     ymax = 5.0*psep
     ymin = 1.0*psep
     ymax = 2.0*psep
     labely = 'h'
  endif
!
!--set up y axis
!
  dy = (ymax - ymin)/REAL(ny-1)
  do i=1,ny
     yaxis(i) = (ymin + (i-1)*dy)
  enddo 
!
!--set up an array of wavenumber values
!
  kxmin = 0.0
  kxmax = 2.*pi
  dkx = (kxmax - kxmin)/REAL(nkx)
  do i=1,nkx
     kxarray(i) = kxmin + i*dkx
  enddo

  mainloop: do iloop=1,nplots
     if (nplots.gt.1) iplotpos = iloop
!     neps = neps + 1
     print*,'eps = ',eps,' neps = ',neps
     write(title,"(a,f3.1,a,i1)") '\ge = ',eps,', n = ',neps
!
!--calculates the 1D dispersion relation for SPH (no av)
!  
  do j=1,ny
     if (negstress) then  ! y axis is R
	R = yaxis(j)
     else                 ! y axis is h
        h = yaxis(j)
     endif
     do i=1,nkx
        kx = kxarray(i)/psep
        omegasq = omegasq1D(h,kxarray(i),cs2,gamma,pmass,rhozero,R,eps,neps,x,npart,ipart)
        dat(i,j) = omegasq/(kx**2*cs2)
!!       print*,i,j,' kx = ',kx,' h = ',h,' R = ',R,' dat = ',dat(i,j)
        if (i.eq.ikx) then
	   tterm1(j) = term1/(kx**2*cs2)
	   tterm2(j) = term2/(kx**2*cs2)
	   tterm3(j) = term3/(kx**2*cs2)
	   tterm4(j) = term4/(kx**2*cs2)
           rhosum(j) = rhosumi
	endif
     enddo
  enddo
!
!--now plot this using PGPLOT
!
  if (contours) then
!!  call pgbegin(0,'?',1,1)
     if (nacross*ndown.eq.1) call pgsch(1.2)
     if (mod(iplotpos,nacross*ndown).eq.1 .or. nacross*ndown.eq.1) call pgpage

     call setpage2(iplotpos,nacross,ndown,kxmin,kxmax,ymin,ymax-0.001,'kx',trim(labely), &
                   ' ',0,0,&
                   0.,0.,0.,0.,0.,0.,.false.,.true.)

     !call danpgtile(iplotpos,nacross,ndown,kxmin,kxmax,ymin,ymax-0.001, &  ! tiled plots
     !            'kx',TRIM(labely),' ',0,0)
!
! plot eps/neps/other labels on plot
!
     write(text,"(a,f3.1)") '\ge = ',eps
!     call pgmtxt('T',-1.5,0.96,1.0,trim(text))
     write(text,"(a,i1)") 'n = ',neps     
!     call pgmtxt('T',-3.0,0.96,1.0,trim(text))
!!     call pgmtxt('T',-3.0,0.98,1.0,trim(kernelname))
!!     call pgmtxt('T',-1.5,0.98,1.0,'Morris'' formalism')

!  call pgenv(kxmin,kxmax,ymin,ymax,0,0)            ! use this for movie
!  call pglabel('kx',TRIM(labely),TRIM(kernelname)) ! use this for movie

     trans(1) = kxmin - 0.5*dkx               ! this is for the pgimag call
     trans(2) = dkx                  ! see help for pgimag/pggray/pgcont
     trans(3) = 0.0
     trans(4) = ymin - 0.5*dy
     trans(5) = 0.0
     trans(6) = dy
!
!--set contour levels
! 
     datmin = minval(dat)
     datmax = maxval(dat)
     print*,'datmax = ',datmax,datmin  !!maxval(dat)*100,nint(maxval(dat)*100)
     dcont = (datmax-datmin)/real(ncont)   ! even contour levels
     do i=1,ncont
        levels(i) = datmin + real(i)*dcont
     enddo
     call pgcont(dat,nkx,ny,1,nkx,1,ny,levels,ncont,trans)
!
!--label levels
!
     do i=1,ncont
        write(string,"(f5.2)") levels(i)
        !!print*,'level = ',levels(i),string
        call pgqch(charheight) ! query character height
        call pgsch(0.5*charheight)   ! shrink character height
	!!call pgsch(0.5)
        call pgconl(dat,nkx,ny,1,nkx,1,ny,levels(i),trans,TRIM(string),35,14)
        call pgsch(charheight) ! restore character height
     enddo
  
  else
!
!--plot just the line along kx = 0
!
     !!call pgsch(1.2)
     call pgsch(1.0)
     !!if (iplotpos.eq.1) call pgpage
     datmin = 0.5  !!min(minval(tterm1),minval(tterm2)) !!,minval(tterm3),minval(tterm4)) !-2.0   !!1.0
     datmax = 1.5  !!max(maxval(tterm1),maxval(tterm2))   !,maxval(tterm3),maxval(tterm4))
     !datmax = maxval(dat(ikx,1:ny))   !!5.5 !!sqrt(maxval(dat(ikx,1:ny)))

     call setpage2(iplotpos,nacross,ndown,ymin,ymax,datmin,datmax,labely,'norm', &
                   ' ',0,0,&
                   0.,0.,0.,0.,0.,0.,.false.,.true.)
!     call danpgtile(iplotpos,nacross,ndown, &
!          ymin,ymax,datmin,datmax,labely,'cs',' ',0,0)
     if (nplots.gt.1) call pgsls(iplotpos)
!     call pgline(ny,yaxis(1:ny),sqrt(abs(dat(ikx,1:ny))))
     !call pgline(ny,yaxis(1:ny),dat(ikx,1:ny))
!
!--plot contributions from first and second derivatives
!
     !--set legend position
     ipos = 1
     hpos = 0.1   ! horizontal position as % of viewport
     vpos = 1.5  ! vertical position in character heights from top
     call pgsls(2) !sci(2)
     call pgsci(2)
!--second derivative normalisation term
     call pgline(ny,yaxis(1:ny),0.5*tterm1(1:ny))
     print*,' tterm = ',0.5*tterm1(1:10),ny
     if (iplotpos.eq.1) call legend(ipos,'grad^2 W norm',hpos,vpos)
     call pgsls(3)
     call pgline(ny,yaxis(1:ny),tterm2(1:ny)+2.0)
     if (iplotpos.eq.1) call legend(ipos,'grad W norm',hpos,vpos)
!     call pgsls(4)
!     call pgline(ny,yaxis(1:ny),tterm3(1:ny))
!     if (iplotpos.eq.1) call legend(ipos,'aniso dW',hpos,vpos)
!     call pgsls(5)
!     call pgline(ny,yaxis(1:ny),tterm4(1:ny))
!     if (iplotpos.eq.1) call legend(ipos,'aniso d\u2\dW',hpos,vpos)
!     call pgsci(1)
!     call pgsls(1)
!
!--also plot density estimate
!
     call pgsls(1)
     call pgsci(2)
     call pgline(ny,yaxis(1:ny),rhosum(1:ny))
     if (iplotpos.eq.1) call legend(ipos,'density estimate',hpos,vpos)
     call pgsci(1)
     read*

     write(text,"(a,f3.1,a,i1)") '\ge = ',eps,', n = ',neps
     !if (nplots.gt.1) then
     !   call legend(iplot,trim(text),hpos,vpos)   
     !else
     !   call legend(0,'cubic spline, h=1.2\gDp',hpos,vpos)  
     !endif
     !
     !--plot line from file
     !
!     if (iplot.eq.1) then
!        open(unit=11,file='legend',status='old')
!20      continue
!        do i=1,4
!           call getarg(i,filename)
!	   call pgsls(i)
!           if (filename(1:1).ne.' ') then
!	      call plotfreq(filename)
!	      read(11,"(a)",end=21) text
!21            continue
!	      print*,trim(text),': ',filename
!              call legend(i,trim(text),hpos,vpos)
!	   endif
!	enddo
!        close(unit=11)
!     else
!        !!call legend(5,trim(kernelname)//', fixed h',hpos,vpos)     
!     endif

     call pgsch(charheight)
     call pgsls(1)
  endif

  enddo mainloop

end subroutine kernelstability1D
!
! function to calculate the dispersion relation for a given
! smoothing length, wavenumber, sound speed, density and particle mass
! 
! calculates this using the analytic dispersion relation but evaluates
! the summations by summing over an array of particles
! In principle could be used for non-uniform particle setups but I have
! never done so.
!
! arguments:
!           x(npart) : array of particle positions
!           ipart    : position of central particle in the array
!           npart    : number of particles
!           RR       : negative stress parameter (RR=1 gives hydro)
!
real function omegasq1D(hh,kx,cs2,gamma,pmass,rhozero,RR,eps,neps,x,npart,ipart)
  use kernels, only:wij,grwij,grgrwij,grwijalt,grgrwijalt, &
                    ikern,dq2table,ddq2table
  implicit none
  integer :: i,index,index1, ipart, npart, neps
  real, dimension(npart) :: x
  real :: hh,kx,cs2,pmass,gamma,rhozero,RR
  real :: przero,Bxzero,priso,praniso
  real :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8
  real :: dxx, dgrwdx, dgrgrwdx, gradW, gradgradW
  real :: Wjoe, W, dWdx, gradWcorr, gradgradWcorr
  real :: dx, dr, q2, eps, term, term1,term2,term3,term4,rhosumi
  common /terms/ term1,term2,term3,term4,rhosumi

  sum1 = 0.
  sum2 = 0.
  sum3 = 0.
  sum4 = 0.
  sum5 = 0.
  sum6 = 0.
  sum7 = 0.
  sum8 = 0.
  rhosumi = 0.
!
!--calculate pressure from cs2 and gamma
!  
  przero = rhozero*cs2/gamma
!
!--Bx from the input negative stress parameter
!
  Bxzero = sqrt(2.*rhozero*cs2*(1.-RR))
!
!--calculate kernel in denominator of anticlumping term
!
  q2 = (0.2)**2
  index = INT(q2*ddq2table)
  index1 = index + 1
  IF (index.GT.ikern) index = ikern
  IF (index1.GT.ikern) index1 = ikern
  dxx = q2 - index*dq2table  
  dWdx = (wij(index1)-wij(index))*ddq2table
  Wjoe = (wij(index) + dWdx*dxx)/hh 
!!  print*,'wjoe = ',wjoe,q2,index,hh,1./hh
  if (wjoe.lt.1.e-5) then
     wjoe = 1.0
     eps = 0.0
  endif
!
!--calculate summation terms
!
  do i=1,npart
     dx = x(i) - x(ipart)
     if (i.ne.ipart) then
        dr = dx/abs(dx)  ! unit vector in x direction
     else
        dr = 0.
     endif
!
!--find nearest index in kernel table
! 
     q2 = (dx/hh)**2
     index = INT(q2*ddq2table)
     index1 = index + 1
     IF (index.GT.ikern) index = ikern
     IF (index1.GT.ikern) index1 = ikern
     dxx = q2 - index*dq2table
     dgrwdx = (grwij(index1)-grwij(index))*ddq2table
     dgrgrwdx = (grgrwij(index1)-grgrwij(index))*ddq2table
     gradW = (grwij(index) + dgrwdx*dxx)/hh**2
     gradgradW = (grgrwij(index) + dgrgrwdx*dxx)/hh**3 
     
     dWdx = (wij(index1)-wij(index))*ddq2table
     W = (wij(index) + dWdx*dxx)/hh 
     
     !--calculate corrected kernels for anticlumping
     if (neps.ne.0) then
!        term = (1.-0.5*eps*(w/wjoe)**neps)
!        gradWcorr = gradW*term
!        gradgradWcorr = gradgradW*term - 0.5*eps*neps*W**(neps-1)/wjoe**neps*gradW*gradW
         dgrwdx = (grwijalt(index1)-grwijalt(index))*ddq2table
         dgrgrwdx = (grgrwijalt(index1)-grgrwijalt(index))*ddq2table
         gradWcorr = (grwijalt(index) + dgrwdx*dxx)/hh**2
         gradgradWcorr = (grgrwijalt(index) + dgrgrwdx*dxx)/hh**3
     else
        gradWcorr = gradW
	gradgradWcorr = gradgradW
     endif 
     sum1 = sum1 + (1. - COS(kx*dx))*gradgradW
     sum2 = sum2 + SIN(kx*dx)*gradW*dr ! times unit vector in x dir
     sum3 = sum3 + (1. - COS(kx*dx))*gradgradWcorr
     sum4 = sum4 + SIN(kx*dx)*gradWcorr*dr
     sum5 = sum5 + dx*gradW*dr
!     sum6 = sum6 + dx*gradWcorr*dr
     sum7 = sum7 + 0.5*dx**2*gradgradW
!     sum8 = sum8 + 0.5*dx**2*gradgradWcorr
     rhosumi = rhosumi + W
  enddo
! 
!!  if (kx.lt.0.1) print*,sum1,kx**2,sum2/kx,sum3,sum4,sum5,sum6,sum7
!  sum1 = kx**2
!  sum3 = kx**2
!  sum4 = -kx
!  sum2 = sum2/sum5
!  sum4 = sum4/sum6
!  sum1 = sum1/sum7
!  sum3 = sum3/sum8

  priso = przero + 0.5*Bxzero**2
  praniso = -Bxzero**2
!
!--dispersion relation (adiabatic and isothermal)
!  isotropic terms using normal kernel
!
  term1 = (2.*pmass/rhozero)*(priso/rhozero)*sum1
  term2 =  - (pmass/rhozero)**2*(2.*priso/rhozero - cs2)*sum2**2
  
  omegasq1D = term1 + term2
!
!--add anisotropic terms (using anticlumping kernel)
!
  term3 = (2.*pmass/rhozero)*praniso/rhozero*sum3
  term4 = - (pmass/rhozero)**2*(2.*praniso/rhozero)*sum2*sum4

  omegasq1D = omegasq1D + term3 + term4
  
!  if (hh.lt.) then
!     print*,hh,kx,' terms = ',(2.*pmass/rhozero)*(priso/rhozero)*sum1, &
 !            - (pmass/rhozero)**2*(2.*priso/rhozero - cs2)*sum2**2, & 
!              + (2.*pmass/rhozero)*praniso/rhozero*sum3, &
!	      - (pmass/rhozero)**2*(2.*praniso/rhozero)*sum2*sum4,omegasq1D
!  endif
end function omegasq1D

!
! subroutine to plot to content of the .freq file
!
subroutine plotfreq(filename)
 implicit none
 integer, parameter :: max = 100
 integer :: i,npts
 real, dimension(max) :: hfact,freq
 character(len=50) :: filename
 
 if (index(filename,'.freq').eq.0) then
    filename = trim(filename)//'.freq'
 endif
 
 print*,'opening ',trim(filename)
 
 open(1,file=filename,status='old',err=20)
 do i=1,max
    read(1,*,end=10) hfact(i),freq(i)
    !!print*,i,hfact(i),freq(i)
 enddo
10 continue
 close(1)
 npts = i-1

! if (filename(1:3).eq.'hfi') then
    call pgpt(npts,hfact(1:npts),freq(1:npts),17)
! else
!    call pgline(npts,hfact(1:npts),freq(1:npts))
! endif
! call pgpt(npts,hfact(1:npts),freq(1:npts),17)
20 continue

 return
end subroutine plotfreq
