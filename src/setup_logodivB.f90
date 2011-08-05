!!------------------------------------------------------------------------!!
!! Setup for a div B advection problem                                    !!
!!  with a logo as the initial conditions                                 !!
!! reads from a .xpm file as output by GIMP
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer, parameter :: maxlines = 1000
 integer :: i,j,ntot,npartx,nparty,ipart
 integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
 integer :: npixx,npixy,iline
 real :: denszero,przero
 real :: pri,rbump,rr
 real :: totmass,gam1,massp,const
 real :: term,q2,wab,grkern,hfacwab,rr2,hh2,xi,yi,dxpix
 real, dimension(ndim) :: xorigin, dx
 real, dimension(ndimV) :: Bzero
 real, dimension(:,:), allocatable :: dat
 character (len=10) :: dummy
 character (len=1000) :: line
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for this problem'
 if (ndimV.ne.3) stop ' need ndimV=3 for this problem'

 const = sqrt(4.*pi)
!
!--open image file (.xpm)
!
 print*,' opening image file : logo.xpm'
 open(unit=55,file='logo.xpm',status='old',form='formatted')
 !
 !--read header
 !
 read(55,*) dummy
 read(55,*) dummy
 read(55,"(a1,i3,1x,i3)") dummy(1:1),npixx,npixy
 print*,' IMAGE SIZE = ',npixx,npixy
 do iline=1,11
    read(55,*) dummy
 enddo
 !
 !--allocate memory for dat array
 !
 allocate(dat(npixx,npixy))
 !
 !--read image
 !
 do iline=1,maxlines
    !--read line as a character string
    read(55,"(a)",end=66) line
    npixx = len_trim(line)
    !!print*,'line ',iline,' length = ',npixx
    print*,'line = ',trim(line)
    do j=2,len_trim(line)-2
       select case(line(j:j))
          case ('%')
             dat(j,iline) = 1./const
          case ('#','*','=','-','&','@','+','$')
             dat(j,iline) = 0.5/const
          case default
             dat(j,iline) = 0.
       end select
    enddo
    !!print*,'line= ',dat(:,iline)
    !!read*
 enddo
66 continue
 npixy = iline-1
 print*,'npixy = ',npixy
 close(unit=55)
!
!--set boundaries
!            	    
 ibound = 3	! periodic
 nbpts = 0	! no fixed particles
 xmin(:) = 0.0	! unit square
 xmax(1) = 1.0
 dxpix = xmax(1)/npixx
 xmax(2) = npixy*dxpix
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0	! co-ordinates of the centre of the initial blast
 rbump = 0.125		! radius of the initial bump
 Bzero(:) = 0.
 if (imhd.ne.0) Bzero(3) = 1.0/const	! uniform field in Bz direction
 przero = 6.0		! initial pressure
 denszero = 1.0		! ambient density
 
 gam1 = gamma - 1.

 write(iprint,*) 'Two dimensional div B advection problem (LOGO)'
 write(iprint,10) denszero,rbump,Bzero(3),przero
10 format(/,' density  = ',f10.3,', size of bump = ',f6.3,/, &
            ' Initial Bz   = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass in ambient medium
!
 totmass = denszero*product(xmax(:)-xmin(:))
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
!    dx(:) = x(:,ipart)-xorigin(:) 
!    rr = dot_product(dx,dx)
    Bfield(:,ipart) = Bzero(:)
!    if (rr.le.rbump) then
!       Bfield(1,ipart) = (4096.*rr**4 - 128.*rr**2 + 1.)/const
!    endif  
    pmass(ipart) = massp
    dens(ipart) = denszero
    vel(:,ipart) = 0.
    !!vel(1:ndim,ipart) = 1.0
    pri = przero 
    uu(ipart) = pri/(gam1*denszero)
    hh(ipart) = hfact*(pmass(ipart)/dens(ipart))**dndim
 enddo

 print*,'interpolating to particles ...'
!
!--interpolate from pixels to particles
!
 do i=1,npart
    !
    !--find nearest pixel
    !
    ipix = int((x(1,i) - xmin(1))/dxpix)
    jpix = int((xmax(2) - x(2,i))/dxpix)
    Bfield(1,i) = dat(ipix,jpix)
    !
    !--find contributing pixel
    !
!    ipixmin = min(int((x(1,i) - 2.*hh(i) - xmin(1))/dxpix),1)
!    ipixmax = max(int((x(1,i) + 2.*hh(i) - xmin(1))/dxpix),npixx)
!    jpixmin = min(int((x(2,i) - 2.*hh(i) - xmin(2))/dxpix),1)
!    jpixmax = max(int((x(2,i) + 2.*hh(i) - xmin(2))/dxpix),npixy)
    !
    !--interpolate from these pixels to current particle
    !
!    term = pmass(i)/dens(i)
!    hh2 = hh(i)**2
!    hfacwab = 1./hh(i)**2
!    print*,'particle ',i,x(:,i),' pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

!    do jpix=jpixmin,jpixmax
!       yi = xmin(2) + (jpix-1)*dxpix
!       do ipix=ipixmin,ipixmax
!          xi = xmin(1) + (ipix-1)*dxpix
!          rr2 = xi**2 + yi**2
!          q2 = rr2/hh2
!          if (q2.lt.4.) then
!             call interpolate_kernel(q2,wab,grkern)
!             wab = wab/hfacwab
!             
!             Bfield(1,i) = Bfield(1,i) + term*dat(ipix,jpix)*wab
!          endif
!       enddo
!    enddo
 enddo

 deallocate(dat)


!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  Exiting subroutine setup'
            
 return
end
