!----------------------------------------------------------------------
!
!  Module containing all of the routines required for 3D projections
!  (where rendered quantity is integrated along the line of sight)
!
!----------------------------------------------------------------------

module projections3D
 implicit none

 integer, parameter :: maxcoltable = 1000
 real, parameter, private :: dpi = 1./3.1415926536
 real, parameter :: radkernel = 2.0
 real, parameter :: radkernel2 = radkernel*radkernel
 real, dimension(0:maxcoltable) :: coltable
 real, parameter, private :: dq2table = radkernel*radkernel/real(maxcoltable)
 real, parameter, private :: ddq2table = 1./dq2table

 public :: setup_integratedkernel
 public :: interpolate3D_projection
 public :: interpolate3D_proj_vec
 public :: wfromtable

contains

subroutine setup_integratedkernel
!-------------------------------------------------------------
!     tabulates the integral through the cubic spline kernel
!     tabulated in (r/h)**2 so that sqrt is not necessary
!-------------------------------------------------------------
 implicit none
 integer :: i,j
 real :: rxy2,deltaz,dz,z,q2,q,wkern,coldens
 integer, parameter :: npts = 100

 print "(1x,a)",'setting up integrated kernel table...'

 do i=0,maxcoltable
!
!--tabulate for (cylindrical) r**2 between 0 and radkernel**2
!
    rxy2 = i*dq2table
!
!--integrate z between 0 and sqrt(radkernel^2 - rxy^2)
!
    deltaz = sqrt(radkernel2 - rxy2)
    dz = deltaz/real(npts-1)
    coldens = 0.
    
    do j=1,npts
       z = (j-1)*dz
       q2 = rxy2 + z*z
       q = sqrt(q2)
       if (q.le.1.0) then
          wkern = 1. - 1.5*q2 + 0.75*q2*q
       elseif(q.lt.2.0) then
          wkern = 0.25*(2.-q)**3
       else
          wkern = 0.
       endif
       if (j.eq.1 .or. j.eq.npts) then
          coldens = coldens + 0.5*wkern*dz ! trapezoidal rule
       else
          coldens = coldens + wkern*dz
       endif
    enddo
    coltable(i)=2.0*coldens*dpi
 end do
 
 return
end subroutine setup_integratedkernel
!
! This function interpolates from the table of integrated kernel values
! to give w(q)
!
real function wfromtable(q2)
 implicit none
 real :: q2,dxx,dwdx
 integer :: index, index1
 !
 !--find nearest index in table
 !
 index = int(q2*ddq2table)
 index1 = min(index + 1,maxcoltable)
 !
 !--find increment along from this index
 !
 dxx = q2 - index*dq2table
 !
 !--find gradient
 !
 dwdx = (coltable(index1) - coltable(index))*ddq2table
 !
 !--compute value of integrated kernel
 !
 wfromtable = coltable(index) + dwdx*dxx

end function wfromtable

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
! 
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     The (dimensionless) weight for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     the interface is written in this form to avoid floating exceptions
!     on physically scaled data.
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!     3D perspective added Nov 2005
!--------------------------------------------------------------------------

subroutine interpolate3D_projection(x,y,z,hh,weight,dat,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidth,zobserver,dscreen)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  real, intent(in) :: xmin,ymin,pixwidth,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: datsmooth

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: iprintinterval, iprintnext, iprogress, itmin
  real :: hi,hi1,hi21,radkern,wab,rab2,q2
  real :: term,dx,dy,dy2,xpix,ypix,zfrac
  real :: t_start,t_end,t_used,tsec
  logical :: iprintprogress

  datsmooth = 0.
  term = 0.
  print "(1x,a)",'projecting from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate3D_proj: error: pixel width <= 0'
     return
  endif
  !
  !--check column density table has actually been setup
  !
  if (abs(coltable(1)).le.1.e-5) then
     call setup_integratedkernel
  endif
  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
!
!--get starting CPU time
!
  call cpu_time(t_start)
  
  over_particles: do i=1,npart
     !
     !--report on progress
     !
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     term = weight(i)*hi*dat(i) ! h gives the z length scale (NB: no perspective)
     if (hi.le.0.) then
        print*,'interpolate3D_proj: error: h <= 0 ',i,hi
        return
     elseif (abs(dscreen).gt.tiny(dscreen)) then
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif
     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = 2.*hi  !radius of the smoothing kernel
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidth) + 1

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
     if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix-0.5)*pixwidth
        dy = ypix - y(i)
        dy2 = dy*dy
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix-0.5)*pixwidth
           dx = xpix - x(i)
           rab2 = dx*dx + dy2
           q2 = rab2*hi21
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !                  
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab          
           endif

        enddo
     enddo

  enddo over_particles
!
!--get ending CPU time
!
  call cpu_time(t_end)
  t_used = t_end - t_start
  if (t_used.gt.60) then
     itmin = int(t_used/60)
     tsec = t_used - (itmin*60)
     print*,'completed in ',itmin,' min ',tsec,'s'
  else
     print*,'completed in ',t_used,'s'
  endif
  
  return

end subroutine interpolate3D_projection

!--------------------------------------------------------------------------
!
!     Same as previous but for a vector quantity
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price 23/12/04
!--------------------------------------------------------------------------

subroutine interpolate3D_proj_vec(x,y,z,hh,weight,vecx,vecy,npart,&
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidth,zobserver,dscreen)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  real, intent(in) :: xmin,ymin,pixwidth,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,hi21,radkern,q2,wab,rab2,const,zfrac
  real :: termx,termy,dx,dy,dy2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  print "(1x,a)",'projecting vector from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_proj_vec: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !      
  over_particles: do i=1,npart
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     const = weight(i)*hi ! h gives the z length scale (NB: no perspective)
     if (hi.le.0.) then
        print*,'interpolate3D_proj_vec: error: h <= 0 ',i,hi
        return
     elseif (abs(dscreen).gt.tiny(dscreen)) then
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif
     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = 2.*hi    ! radius of the smoothing kernel
        
     termx = const*vecx(i)
     termy = const*vecy(i)
     !
     !--for each particle work out which pixels it contributes to
     !               
     ipixmin = int((x(i) - radkern - xmin)/pixwidth)
     jpixmin = int((y(i) - radkern - ymin)/pixwidth)
     ipixmax = int((x(i) + radkern - xmin)/pixwidth) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidth) + 1

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix-0.5)*pixwidth
        dy = ypix - y(i)
        dy2 = dy*dy
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix-0.5)*pixwidth
           dx = xpix - x(i)
           rab2 = dx**2 + dy2
           q2 = rab2*hi21
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !  
              vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
              vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
           endif

        enddo
     enddo

  enddo over_particles

  return

end subroutine interpolate3D_proj_vec

end module projections3D
